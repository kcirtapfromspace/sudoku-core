use super::*;

#[test]
fn test_box_index() {
    assert_eq!(Position::new(0, 0).box_index(), 0);
    assert_eq!(Position::new(0, 3).box_index(), 1);
    assert_eq!(Position::new(3, 0).box_index(), 3);
    assert_eq!(Position::new(8, 8).box_index(), 8);
    assert_eq!(Position::new(4, 4).box_index(), 4);
}

#[test]
fn test_diagonals() {
    assert!(Position::new(0, 0).is_on_main_diagonal(9));
    assert!(Position::new(4, 4).is_on_main_diagonal(9));
    assert!(Position::new(8, 8).is_on_main_diagonal(9));
    assert!(!Position::new(0, 1).is_on_main_diagonal(9));

    assert!(Position::new(0, 8).is_on_anti_diagonal(9));
    assert!(Position::new(4, 4).is_on_anti_diagonal(9));
    assert!(Position::new(8, 0).is_on_anti_diagonal(9));
    assert!(!Position::new(0, 0).is_on_anti_diagonal(9));
}
