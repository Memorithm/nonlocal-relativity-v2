use scirust_nonlocal_relativity::{
    CylindricalMinkowski, DiscreteConnectionTransport, HistoryEntry, HistoryTransport,
    WorldlineState, transport_vector_along_polyline,
};
use scirust_relativity::{Connection, Minkowski};

fn identity<const D: usize>() -> [[f64; D]; D] {
    let mut matrix = [[0.0_f64; D]; D];
    for i in 0..D {
        matrix[i][i] = 1.0;
    }
    matrix
}

fn matrix_vector<const D: usize>(matrix: &[[f64; D]; D], vector: &[f64; D]) -> [f64; D] {
    let mut out = [0.0_f64; D];
    for row in 0..D {
        for column in 0..D {
            out[row] += matrix[row][column] * vector[column];
        }
    }
    out
}

fn matrix_multiply<const D: usize>(left: &[[f64; D]; D], right: &[[f64; D]; D]) -> [[f64; D]; D] {
    let mut out = [[0.0_f64; D]; D];
    for row in 0..D {
        for column in 0..D {
            for inner in 0..D {
                out[row][column] += left[row][inner] * right[inner][column];
            }
        }
    }
    out
}

/// Build the exact linear operator of the *existing* segment transport by
/// probing it on the canonical basis. This intentionally does not optimize
/// production code yet: it is an independent characterization oracle proving
/// that the existing Heun segment transport is a linear map that can be
/// factored and composed before the optimized implementation lands.
fn probed_segment_operator<B, const D: usize>(
    background: &B,
    from: &WorldlineState<D>,
    to: &WorldlineState<D>,
    step: f64,
) -> [[f64; D]; D]
where
    B: Connection<D>,
{
    let transport = DiscreteConnectionTransport;
    let mut operator = [[0.0_f64; D]; D];

    for column in 0..D {
        let mut basis = [0.0_f64; D];
        basis[column] = 1.0;
        let transported = transport
            .transport_segment(column, background, basis, from, to, step)
            .expect("basis transport must succeed");
        for row in 0..D {
            operator[row][column] = transported[row];
        }
    }

    operator
}

fn max_abs_difference<const D: usize>(left: &[f64; D], right: &[f64; D]) -> f64 {
    let mut max = 0.0_f64;
    for component in 0..D {
        max = max.max((left[component] - right[component]).abs());
    }
    max
}

#[test]
fn zero_connection_operator_is_bit_exact_identity() {
    let from = WorldlineState::new([0.0, 1.0, 2.0, 3.0], [1.2, 0.2, -0.1, 0.05]);
    let to = WorldlineState::new([0.05, 1.01, 1.99, 3.02], [1.18, 0.21, -0.08, 0.04]);
    let operator = probed_segment_operator(&Minkowski, &from, &to, 0.05);

    assert_eq!(operator, identity::<4>());
}

#[test]
fn segment_transport_is_reconstructed_by_one_linear_operator() {
    let background = CylindricalMinkowski;
    let transport = DiscreteConnectionTransport;
    let from = WorldlineState::new([0.0, 5.0, 0.7, 0.0], [1.2, 0.15, 0.08, -0.05]);
    let to = WorldlineState::new([0.03, 5.01, 0.72, -0.001], [1.19, 0.1515, 0.079, -0.05]);
    let step = 0.03;
    let operator = probed_segment_operator(&background, &from, &to, step);

    let vectors = [
        [1.3, 0.2, -0.1, 0.05],
        [0.7, -0.4, 0.3, 0.2],
        [-0.2, 0.6, 0.1, -0.8],
        [2.0, 0.0, 0.0, 0.0],
    ];

    for (index, vector) in vectors.into_iter().enumerate() {
        let direct = transport
            .transport_segment(index, &background, vector, &from, &to, step)
            .expect("direct segment transport must succeed");
        let factored = matrix_vector(&operator, &vector);
        let error = max_abs_difference(&direct, &factored);
        assert!(
            error <= 2.0e-15,
            "operator reconstruction error {error:.3e} for vector {index}: direct={direct:?}, factored={factored:?}"
        );
    }
}

#[test]
fn composed_segment_operators_reproduce_polyline_transport() {
    let background = CylindricalMinkowski;
    let waypoints = [
        HistoryEntry::new([0.00, 5.000, 0.700, 0.000], [1.20, 0.150, 0.080, -0.050], 0.00),
        HistoryEntry::new([0.03, 5.010, 0.720, -0.001], [1.19, 0.151, 0.079, -0.050], 0.03),
        HistoryEntry::new([0.07, 5.025, 0.748, -0.003], [1.18, 0.153, 0.078, -0.049], 0.07),
        HistoryEntry::new([0.12, 5.045, 0.785, -0.006], [1.17, 0.154, 0.077, -0.048], 0.12),
    ];

    let mut prefix = identity::<4>();
    for window in waypoints.windows(2) {
        let from = WorldlineState::new(window[0].coordinates, window[0].velocity);
        let to = WorldlineState::new(window[1].coordinates, window[1].velocity);
        let step = window[1].parameter - window[0].parameter;
        let segment = probed_segment_operator(&background, &from, &to, step);
        prefix = matrix_multiply(&segment, &prefix);
    }

    let vectors = [
        [1.3, 0.2, -0.1, 0.05],
        [0.7, -0.4, 0.3, 0.2],
        [-0.2, 0.6, 0.1, -0.8],
    ];

    for vector in vectors {
        let direct = transport_vector_along_polyline(
            &background,
            &DiscreteConnectionTransport,
            vector,
            &waypoints,
        )
        .expect("polyline transport must succeed");
        let factored = matrix_vector(&prefix, &vector);
        let error = max_abs_difference(&direct, &factored);
        assert!(
            error <= 4.0e-15,
            "prefix factorization error {error:.3e}: direct={direct:?}, factored={factored:?}"
        );
    }
}

#[test]
fn operator_construction_and_composition_are_deterministic_bit_for_bit() {
    let background = CylindricalMinkowski;
    let from = WorldlineState::new([0.0, 5.0, 0.7, 0.0], [1.2, 0.15, 0.08, -0.05]);
    let to = WorldlineState::new([0.03, 5.01, 0.72, -0.001], [1.19, 0.1515, 0.079, -0.05]);

    let first = probed_segment_operator(&background, &from, &to, 0.03);
    let second = probed_segment_operator(&background, &from, &to, 0.03);
    assert_eq!(first, second);

    let first_square = matrix_multiply(&first, &first);
    let second_square = matrix_multiply(&second, &second);
    assert_eq!(first_square, second_square);
}
