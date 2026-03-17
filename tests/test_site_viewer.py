from gui.ui.widgets.site_viewer import SiteViewer


def test_site_viewer_displays_expected_cells(qt_app):
    records = [
        ("spA", "AAAA"),
        ("spB", "AAAA"),
        ("spC", "BBBB"),
        ("spD", "BBBB"),
        ("spE", "BBBB"),
        ("spF", "CCCC"),
    ]

    viewer = SiteViewer(
        records,
        ["spA", "spB"],
        ["spC", "spD"],
        ["spE"],
        gene_name="demo",
    )
    viewer.show()
    qt_app.processEvents()

    top_model = viewer.top_right_table.model()
    bottom_model = viewer.bottom_right_table.model()

    assert viewer.top_right_table.rowCount() == 13
    assert viewer.top_right_table.columnCount() == 4
    assert top_model.index(0, 0).data() == "1"
    assert top_model.index(1, 0).data() == "4*"
    assert top_model.index(4, 0).data() == "A"
    assert top_model.index(8, 0).data() == "B"
    assert top_model.index(12, 0).data() == "B"

    assert viewer.bottom_right_table.rowCount() == 2
    assert viewer.bottom_right_table.columnCount() == 4
    assert bottom_model.index(1, 0).data() == "C"

    viewer.close()
