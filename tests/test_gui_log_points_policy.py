from gui.core.config import ESLConfig
from gui.ui.pages.parameters_page import ParametersPage


def test_hidden_mode_auto_points_follow_output_selection(qt_app):
    cfg = ESLConfig()
    params = ParametersPage(cfg)
    params.show()

    assert not params.show_advanced_chk.isChecked()
    assert params.num_log_points.value() == 20

    params.genes_only_btn.setChecked(True)
    assert params.num_log_points.value() == 4
    assert cfg.num_points == 4

    params.preds_only_btn.setChecked(True)
    assert params.num_log_points.value() == 20
    assert cfg.num_points == 20

    params.close()


def test_manual_points_stick_after_hiding_advanced(qt_app):
    cfg = ESLConfig()
    params = ParametersPage(cfg)
    params.show()

    params.show_advanced_chk.setChecked(True)
    params.num_log_points.setValue(37)
    assert cfg.num_points == 37

    params.show_advanced_chk.setChecked(False)
    assert params.num_log_points.value() == 37

    params.genes_only_btn.setChecked(True)
    assert params.num_log_points.value() == 37

    params.both_outputs_btn.setChecked(True)
    assert params.num_log_points.value() == 37

    params.close()


def test_hiding_advanced_without_manual_points_edit_reapplies_policy(qt_app):
    cfg = ESLConfig()
    params = ParametersPage(cfg)
    params.show()

    params.genes_only_btn.setChecked(True)
    assert params.num_log_points.value() == 4

    params.show_advanced_chk.setChecked(True)
    params.both_outputs_btn.setChecked(True)
    # Because the user did not manually edit the field, output-mode changes should
    # still keep the visible value in sync even while advanced options are open.
    assert params.num_log_points.value() == 20

    params.show_advanced_chk.setChecked(False)
    assert params.num_log_points.value() == 20
    assert cfg.num_points == 20

    params.close()


def test_advanced_visible_manual_override_blocks_auto_switch(qt_app):
    cfg = ESLConfig()
    params = ParametersPage(cfg)
    params.show()

    params.show_advanced_chk.setChecked(True)
    params.num_log_points.setValue(37)
    assert cfg.num_points == 37

    params.genes_only_btn.setChecked(True)
    assert params.num_log_points.value() == 37
    assert cfg.num_points == 37

    params.close()
