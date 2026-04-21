from PySide6.QtWidgets import QLineEdit

from gui.core.config import ESLConfig
from gui.core.paths import default_output_dir
from gui.ui.main_window import ESLWizard
from gui.ui.pages.parameters_page import ParametersPage
from gui.ui.pages.command_page import CommandPage


def test_checkbox_with_continuous_pheno(qt_app, tmp_path):
    cfg = ESLConfig()
    cfg.output_dir = str(tmp_path / "out")
    cfg.species_pheno_is_continuous = True
    cfg.species_phenotypes_file = "dummy.csv"
    params = ParametersPage(cfg)
    params.update_output_options_state()
    params.show()
    assert params.use_continuous_chk.isVisible()
    assert params.use_continuous_chk.isEnabled()
    params.use_continuous_chk.setChecked(True)
    assert cfg.use_continuous_phenotypes
    assert params.continuous_plot_chk.isVisible()
    params.continuous_plot_chk.setChecked(True)
    assert cfg.make_continuous_plot
    cmd_page = CommandPage(cfg)
    cmd_page.on_enter()
    cmd_page.show()
    qt_app.processEvents()
    assert "--use_continuous_phenotypes" in cmd_page.cmd_display.toPlainText()
    assert "--make_continuous_plot" in cmd_page.cmd_display.toPlainText()
    texts = [w.text() for w in cmd_page.summary_group.findChildren(QLineEdit)]
    assert any("Continuous phenotypes" in t for t in texts)


def test_checkbox_locked_for_response_dir(qt_app, tmp_path):
    cfg = ESLConfig()
    cfg.output_dir = str(tmp_path / "out")
    cfg.response_dir = "respdir"
    cfg.response_matrices_are_continuous = True
    cfg.use_continuous_phenotypes = True
    params = ParametersPage(cfg)
    params.update_output_options_state()
    params.show()
    assert params.use_continuous_chk.isVisible()
    assert not params.use_continuous_chk.isEnabled()
    assert params.use_continuous_chk.isChecked()
    assert not params.continuous_plot_chk.isVisible()
    cmd_page = CommandPage(cfg)
    cmd_page.on_enter()
    cmd_page.show()
    qt_app.processEvents()
    assert "--use_continuous_phenotypes" in cmd_page.cmd_display.toPlainText()
    texts = [w.text() for w in cmd_page.summary_group.findChildren(QLineEdit)]
    assert any("Continuous phenotypes" in t for t in texts)


def test_continuous_toggle_not_hidden_by_advanced_toggle(qt_app):
    cfg = ESLConfig()
    cfg.species_pheno_is_continuous = True
    cfg.species_phenotypes_file = "dummy.csv"
    params = ParametersPage(cfg)
    params.update_output_options_state()
    params.show()
    qt_app.processEvents()

    assert params.use_continuous_chk.isVisible()
    params.show_advanced_chk.setChecked(True)
    qt_app.processEvents()
    assert params.use_continuous_chk.isVisible()
    params.show_advanced_chk.setChecked(False)
    qt_app.processEvents()
    assert params.use_continuous_chk.isVisible()


def test_command_page_prompts_for_output_dir_when_unset(qt_app):
    cfg = ESLConfig()
    cmd_page = CommandPage(cfg)
    cmd_page.on_enter()
    cmd_page.show()
    qt_app.processEvents()

    text = cmd_page.cmd_display.toPlainText()
    assert "--output_dir" in text
    assert default_output_dir() in text


def test_disable_ec_checkbox_updates_command_preview(qt_app, tmp_path):
    cfg = ESLConfig()
    cfg.output_dir = str(tmp_path / "out")
    params = ParametersPage(cfg)
    params.show()
    qt_app.processEvents()

    params.show_advanced_chk.setChecked(True)
    qt_app.processEvents()

    assert params.disable_ec_chk.isVisible()
    assert cfg.disable_ec

    text_default = CommandPage(cfg)
    text_default.on_enter()
    text_default.show()
    qt_app.processEvents()
    assert "--enable_ec" not in text_default.cmd_display.toPlainText()

    params.disable_ec_chk.setChecked(False)
    assert not cfg.disable_ec

    cmd_page = CommandPage(cfg)
    cmd_page.on_enter()
    cmd_page.show()
    qt_app.processEvents()

    text = cmd_page.cmd_display.toPlainText()
    assert "--enable_ec" in text


def test_output_dir_selection_uses_empty_folder_directly(qt_app, tmp_path):
    cfg = ESLConfig()
    params = ParametersPage(cfg)
    empty_dir = tmp_path / "empty_output"
    empty_dir.mkdir()

    resolved, uses_subdir = params._resolve_output_dir_selection(str(empty_dir))

    assert resolved == str(empty_dir)
    assert not uses_subdir


def test_parameters_page_starts_with_default_output_dir(qt_app):
    cfg = ESLConfig()
    params = ParametersPage(cfg)

    assert cfg.output_dir == default_output_dir()
    assert params.output_dir_edit.text() == cfg.output_dir
    assert params.isComplete()


def test_parameters_page_restores_default_output_dir_when_blank(qt_app):
    cfg = ESLConfig()
    cfg.output_dir = ""
    params = ParametersPage(cfg)

    params.update_ui_from_config()

    assert cfg.output_dir == default_output_dir()
    assert params.output_dir_edit.text() == cfg.output_dir
    assert params.isComplete()


def test_wizard_next_enabled_on_parameters_page_with_default_output_dir(qt_app):
    wiz = ESLWizard()
    wiz.setCurrentId(wiz.pageIds()[1])
    qt_app.processEvents()

    next_btn = wiz.button(wiz.WizardButton.NextButton)
    assert wiz.currentPage() is wiz.params_page
    assert next_btn.isEnabled()
    wiz.close()


def test_output_dir_selection_uses_base_name_subdir_for_nonempty_folder(qt_app, tmp_path):
    cfg = ESLConfig()
    params = ParametersPage(cfg)
    parent_dir = tmp_path / "parent_output"
    parent_dir.mkdir()
    (parent_dir / "existing.txt").write_text("x", encoding="utf-8")
    params.output_file_base_name.setText("custom_run")

    resolved, uses_subdir = params._resolve_output_dir_selection(str(parent_dir))

    assert resolved == str(parent_dir / "custom_run")
    assert uses_subdir


def test_output_dir_selection_tracks_base_name_for_nonempty_parent(qt_app, tmp_path):
    cfg = ESLConfig()
    params = ParametersPage(cfg)
    parent_dir = tmp_path / "parent_output"
    parent_dir.mkdir()
    (parent_dir / "existing.txt").write_text("x", encoding="utf-8")

    params.output_file_base_name.setText("first_run")
    params._apply_output_dir_selection(str(parent_dir))
    assert cfg.output_dir == str(parent_dir / "first_run")

    params.output_file_base_name.setText("renamed_run")
    assert cfg.output_dir == str(parent_dir / "renamed_run")
