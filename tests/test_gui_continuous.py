from PySide6.QtWidgets import QLineEdit

from gui.core.config import ESLConfig
from gui.core.paths import default_output_dir
from gui.ui.main_window import ESLWizard
from gui.ui.pages.parameters_page import ParametersPage
from gui.ui.pages.command_page import CommandPage
from gui.ui.widgets.results_display import PredictionMetricsDialog


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


def test_pheno_names_are_included_in_prediction_command(qt_app, tmp_path):
    cfg = ESLConfig()
    cfg.output_dir = str(tmp_path / "out")
    cfg.pheno_name1 = "C4"
    cfg.pheno_name2 = "C3"
    cfg.no_pred_output = False

    args = cfg.to_cli_args()
    idx = args.index("--pheno_names")
    assert args[idx + 1:idx + 3] == ["C4", "C3"]

    cmd_page = CommandPage(cfg)
    cmd_page.on_enter()
    cmd_page.show()
    qt_app.processEvents()

    text = cmd_page.cmd_display.toPlainText()
    assert "--pheno_names C4 C3" in text


def test_pheno_names_are_omitted_when_predictions_are_disabled():
    cfg = ESLConfig()
    cfg.no_pred_output = True
    cfg.pheno_name1 = "C4"
    cfg.pheno_name2 = "C3"

    assert "--pheno_names" not in cfg.to_cli_args()


def test_prediction_metrics_uses_configured_phenotype_names(qt_app, tmp_path):
    pred_csv = tmp_path / "predictions.csv"
    pred_csv.write_text(
        "species,SPS,num_genes,input_RMSE,true_phenotype\n"
        "Maize,0.8,10,0.1,1\n"
        "Rice,-0.7,10,0.2,-1\n",
        encoding="utf-8",
    )
    cfg = ESLConfig()
    cfg.pheno_name1 = "C4"
    cfg.pheno_name2 = "C3"
    cfg.species_pheno_is_binary = True

    dialog = PredictionMetricsDialog(str(pred_csv), cfg)
    qt_app.processEvents()

    assert dialog._ready
    row_headers = [
        dialog.rows_table.horizontalHeaderItem(i).text()
        for i in range(dialog.rows_table.columnCount())
    ]
    assert "TPR (C4)" in row_headers
    assert "TNR (C3)" in row_headers

    species_headers = [
        dialog.species_table.horizontalHeaderItem(i).text()
        for i in range(dialog.species_table.columnCount())
    ]
    true_pheno_col = species_headers.index("True Phenotype")
    true_pheno_values = {
        dialog.species_table.item(row, true_pheno_col).text()
        for row in range(dialog.species_table.rowCount())
    }
    assert true_pheno_values == {"C4", "C3"}
    dialog.close()


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

    params._apply_output_dir_selection(str(empty_dir))

    assert cfg.output_dir == str(empty_dir)
    assert params.output_dir_edit.text() == str(empty_dir)


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


def test_output_dir_selection_uses_selected_folder_directly_for_nonempty_folder(qt_app, tmp_path):
    cfg = ESLConfig()
    params = ParametersPage(cfg)
    parent_dir = tmp_path / "parent_output"
    parent_dir.mkdir()
    (parent_dir / "existing.txt").write_text("x", encoding="utf-8")
    params.output_file_base_name.setText("custom_run")

    params._apply_output_dir_selection(str(parent_dir))

    assert cfg.output_dir == str(parent_dir)
    assert params.output_dir_edit.text() == str(parent_dir)


def test_output_dir_selection_does_not_track_base_name_changes(qt_app, tmp_path):
    cfg = ESLConfig()
    params = ParametersPage(cfg)
    parent_dir = tmp_path / "parent_output"
    parent_dir.mkdir()
    (parent_dir / "existing.txt").write_text("x", encoding="utf-8")

    params.output_file_base_name.setText("first_run")
    params._apply_output_dir_selection(str(parent_dir))
    assert cfg.output_dir == str(parent_dir)

    params.output_file_base_name.setText("renamed_run")
    assert cfg.output_dir == str(parent_dir)


def test_output_dir_field_is_editable_and_updates_config(qt_app, tmp_path):
    cfg = ESLConfig()
    params = ParametersPage(cfg)
    custom_dir = tmp_path / "typed_output"

    assert not params.output_dir_edit.isReadOnly()

    params.output_dir_edit.setText(str(custom_dir))
    qt_app.processEvents()

    assert cfg.output_dir == str(custom_dir)
