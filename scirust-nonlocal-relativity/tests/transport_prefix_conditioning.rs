use scirust_nonlocal_relativity::{
    CylindricalMinkowski, DiscreteConnectionTransport, HistoryEntry,
    transport_prefix_diagnostics, transport_vector_along_polyline,
};
use scirust_relativity::{Connection, Minkowski};

fn identity<const D: usize>() -> [[f64; D]; D] {
    let mut matrix = [[0.0_f64; D]; D];
    for (index, row) in matrix.iter_mut().enumerate()
    {
        row[index] = 1.0;
    }
    matrix
}

fn matrix_bits_equal<const D: usize>(left: &[[f64; D]; D], right: &[[f64; D]; D]) -> bool {
    left.iter().zip(right.iter()).all(|(left_row, right_row)| {
        left_row
            .iter()
            .zip(right_row.iter())
            .all(|(left_value, right_value)| left_value.to_bits() == right_value.to_bits())
    })
}

fn max_abs_difference<const D: usize>(left: &[f64; D], right: &[f64; D]) -> f64 {
    left.iter()
        .zip(right.iter())
        .map(|(left, right)| (left - right).abs())
        .fold(0.0_f64, f64::max)
}

#[test]
fn minkowski_prefixes_are_bit_exact_identity_with_unit_condition() {
    let waypoints = [
        HistoryEntry::new([0.0, 1.0, 2.0, 3.0], [1.0, 0.1, 0.0, 0.0], 0.0),
        HistoryEntry::new([0.1, 1.01, 2.0, 3.0], [1.0, 0.1, 0.0, 0.0], 0.1),
        HistoryEntry::new([0.2, 1.02, 2.0, 3.0], [1.0, 0.1, 0.0, 0.0], 0.2),
    ];

    let diagnostics =
        transport_prefix_diagnostics(&Minkowski, &waypoints).expect("Minkowski prefixes");
    let expected = identity::<4>();

    assert_eq!(diagnostics.len(), waypoints.len());
    for (segment_count, diagnostic) in diagnostics.iter().enumerate()
    {
        assert_eq!(diagnostic.segment_count, segment_count);
        assert!(matrix_bits_equal(&diagnostic.prefix_operator, &expected));
        assert!(matrix_bits_equal(
            diagnostic.inverse_operator.as_ref().expect("identity inverse"),
            &expected
        ));
        assert_eq!(diagnostic.condition_inf, Some(1.0));
        assert_eq!(diagnostic.min_abs_pivot, Some(1.0));
    }
}

#[test]
fn cylindrical_prefix_matches_scalar_polyline_and_inverts_round_trip() {
    let background = CylindricalMinkowski;
    let waypoints = [
        HistoryEntry::new(
            [0.00, 5.000, 0.700, 0.000],
            [1.20, 0.150, 0.080, -0.050],
            0.00,
        ),
        HistoryEntry::new(
            [0.03, 5.010, 0.720, -0.001],
            [1.19, 0.151, 0.079, -0.050],
            0.03,
        ),
        HistoryEntry::new(
            [0.07, 5.025, 0.748, -0.003],
            [1.18, 0.153, 0.078, -0.049],
            0.07,
        ),
        HistoryEntry::new(
            [0.12, 5.045, 0.785, -0.006],
            [1.17, 0.154, 0.077, -0.048],
            0.12,
        ),
    ];
    let vector = [1.3, 0.2, -0.1, 0.05];

    let diagnostics =
        transport_prefix_diagnostics(&background, &waypoints).expect("cylindrical prefixes");
    let final_prefix = diagnostics.last().expect("final prefix");
    let direct = transport_vector_along_polyline(
        &background,
        &DiscreteConnectionTransport,
        vector,
        &waypoints,
    )
    .expect("scalar polyline transport");
    let factored = DiscreteConnectionTransport::apply_segment_operator(
        &final_prefix.prefix_operator,
        vector,
        0,
    )
    .expect("prefix application");
    let inverse = final_prefix
        .inverse_operator
        .as_ref()
        .expect("regular cylindrical prefix");
    let recovered =
        DiscreteConnectionTransport::apply_segment_operator(inverse, factored, 0)
            .expect("inverse prefix application");

    assert!(max_abs_difference(&direct, &factored) <= 6.0e-15);
    assert!(max_abs_difference(&vector, &recovered) <= 6.0e-15);
    assert!(final_prefix.condition_inf.is_some_and(|value| value >= 1.0));
    assert!(final_prefix.min_abs_pivot.is_some_and(|value| value > 0.0));
}

#[derive(Debug, Clone, Copy)]
struct EndpointSingularConnection;

impl Connection<1> for EndpointSingularConnection {
    fn christoffel(&self, coordinates: &[f64; 1]) -> [[[f64; 1]; 1]; 1] {
        if coordinates[0] < 0.5
        {
            [[[0.0]]]
        }
        else
        {
            [[[2.0]]]
        }
    }
}

#[test]
fn singular_prefix_is_reported_without_manufacturing_an_inverse() {
    let waypoints = [
        HistoryEntry::new([0.0], [1.0], 0.0),
        HistoryEntry::new([1.0], [1.0], 1.0),
    ];

    let diagnostics = transport_prefix_diagnostics(&EndpointSingularConnection, &waypoints)
        .expect("finite singular prefix must still be diagnosable");
    let singular = diagnostics.last().expect("singular prefix");

    assert_eq!(singular.prefix_operator[0][0], 0.0);
    assert!(singular.inverse_operator.is_none());
    assert!(singular.condition_inf.is_none());
    assert_eq!(singular.min_abs_pivot, Some(0.0));
}
