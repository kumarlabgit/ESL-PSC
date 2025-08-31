from PySide6.QtWidgets import QLineEdit

from gui.core.config import ESLConfig
from gui.ui.pages.parameters_page import ParametersPage
from gui.ui.pages.command_page import CommandPage


def test_checkbox_with_continuous_pheno(qt_app):
    cfg = ESLConfig()
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


def test_checkbox_locked_for_response_dir(qt_app):
    cfg = ESLConfig()
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
