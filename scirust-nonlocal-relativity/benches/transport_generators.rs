use criterion::{BenchmarkId, Criterion, black_box, criterion_group, criterion_main};
use scirust_nonlocal_relativity::{
    CylindricalMinkowski, DiscreteConnectionTransport, HistoryTransport, WorldlineState,
};

fn retained_vectors(history_len: usize) -> Vec<[f64; 4]> {
    (0..history_len)
        .map(|index| {
            let scale = index as f64 + 1.0;
            [1.0 + 0.01 * scale, 0.2 / scale, -0.1 / scale, 0.05]
        })
        .collect()
}

fn bench_retained_transport(c: &mut Criterion) {
    let background = CylindricalMinkowski;
    let transport = DiscreteConnectionTransport;
    let from = WorldlineState::new([0.0, 5.0, 0.7, 0.0], [1.2, 0.15, 0.08, -0.05]);
    let to = WorldlineState::new([0.03, 5.01, 0.72, -0.001], [1.19, 0.1515, 0.079, -0.05]);
    let step = 0.03;
    let mut group = c.benchmark_group("retained_transport");

    for history_len in [1_usize, 2, 4, 8, 64, 256]
    {
        let template = retained_vectors(history_len);

        group.bench_with_input(
            BenchmarkId::new("scalar_heun", history_len),
            &history_len,
            |b, _| {
                b.iter(|| {
                    let mut outputs = template.clone();
                    for (index, vector) in outputs.iter_mut().enumerate()
                    {
                        *vector = transport
                            .transport_segment(index, &background, *vector, &from, &to, step)
                            .expect("scalar transport must succeed");
                    }
                    black_box(outputs)
                });
            },
        );

        group.bench_with_input(
            BenchmarkId::new("shared_generator_batch", history_len),
            &history_len,
            |b, _| {
                b.iter(|| {
                    let mut outputs = template.clone();
                    transport
                        .transport_batch(&background, &mut outputs, &from, &to, step)
                        .expect("shared-generator batch transport must succeed");
                    black_box(outputs)
                });
            },
        );
    }

    group.finish();
}

criterion_group!(benches, bench_retained_transport);
criterion_main!(benches);
