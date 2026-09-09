use std::cell::Cell;

use scirust_nonlocal_relativity::{
    CylindricalMinkowski, DiscreteConnectionTransport, HistoryTransport, WorldlineState,
};
use scirust_relativity::Connection;

struct CountingConnection<B> {
    inner: B,
    calls: Cell<usize>,
}

impl<B> CountingConnection<B> {
    const fn new(inner: B) -> Self {
        Self {
            inner,
            calls: Cell::new(0),
        }
    }

    fn calls(&self) -> usize {
        self.calls.get()
    }
}

impl<B, const D: usize> Connection<D> for CountingConnection<B>
where
    B: Connection<D>,
{
    fn christoffel(&self, coordinates: &[f64; D]) -> [[[f64; D]; D]; D] {
        self.calls.set(self.calls.get() + 1);
        self.inner.christoffel(coordinates)
    }
}

fn max_abs_difference<const D: usize>(left: &[f64; D], right: &[f64; D]) -> f64 {
    left.iter()
        .zip(right.iter())
        .fold(0.0_f64, |maximum, (left_value, right_value)| {
            maximum.max((left_value - right_value).abs())
        })
}

#[test]
fn shared_generators_reduce_connection_evaluations_from_two_per_vector_to_two_per_batch() {
    let from = WorldlineState::new([0.0, 5.0, 0.7, 0.0], [1.2, 0.15, 0.08, -0.05]);
    let to = WorldlineState::new([0.03, 5.01, 0.72, -0.001], [1.19, 0.1515, 0.079, -0.05]);
    let step = 0.03;
    let transport = DiscreteConnectionTransport;

    for history_len in [1_usize, 2, 8, 64]
    {
        let template: Vec<[f64; 4]> = (0..history_len)
            .map(|index| {
                let scale = index as f64 + 1.0;
                [1.0 + 0.01 * scale, 0.2 / scale, -0.1 / scale, 0.05]
            })
            .collect();

        let scalar_connection = CountingConnection::new(CylindricalMinkowski);
        let mut scalar_outputs = Vec::with_capacity(history_len);
        for (index, vector) in template.iter().copied().enumerate()
        {
            scalar_outputs.push(
                transport
                    .transport_segment(index, &scalar_connection, vector, &from, &to, step)
                    .expect("scalar Heun transport must succeed"),
            );
        }
        assert_eq!(scalar_connection.calls(), 2 * history_len);

        let batch_connection = CountingConnection::new(CylindricalMinkowski);
        let mut batch_outputs = template.clone();
        transport
            .transport_batch(&batch_connection, &mut batch_outputs, &from, &to, step)
            .expect("shared-generator batch transport must succeed");
        assert_eq!(batch_connection.calls(), 2);

        for (index, (scalar, batch)) in scalar_outputs.iter().zip(batch_outputs.iter()).enumerate()
        {
            let error = max_abs_difference(scalar, batch);
            assert!(
                error <= 3.0e-15,
                "scalar/batch mismatch for history_len={history_len}, vector={index}: {error:.3e}"
            );
        }
    }
}
