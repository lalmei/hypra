import json

from hypra.__main__ import main


def test_backend_info_reports_missing(monkeypatch, capsys):
    monkeypatch.setattr("hypra.__main__.backend_info", lambda **kwargs: {"available": False, "error": "missing"})
    assert main(["backend-info"]) == 2
    assert json.loads(capsys.readouterr().out)["available"] is False


def test_no_expression_eval(tmp_path, capsys):
    request = tmp_path / "request.json"
    request.write_text(json.dumps({"expression": "__import__('os').system('echo unsafe')", "epsilon": {}}))
    assert main(["expand", str(request)]) == 2
    assert json.loads(capsys.readouterr().err)["status"] == "error"


def test_missing_request_file_is_error(tmp_path, capsys):
    assert main(["expand", str(tmp_path / "absent.json")]) == 2
    assert json.loads(capsys.readouterr().err)["status"] == "error"


def test_duplicate_json_keys_rejected(tmp_path, capsys):
    request = tmp_path / "request.json"
    request.write_text('{"through": 1, "through": 2}')
    assert main(["expand", str(request)]) == 2
    assert "Duplicate" in json.loads(capsys.readouterr().err)["message"]


def test_real_typer_application_help_and_backend_info():
    from typer.testing import CliRunner

    from hypra.cli import app

    runner = CliRunner()
    help_result = runner.invoke(app, ['--help'])
    assert help_result.exit_code == 0
    assert 'backend-info' in help_result.stdout
    assert 'expand' in help_result.stdout
    result = runner.invoke(app, ['backend-info'])
    assert result.exit_code == 0
    assert json.loads(result.stdout)['available'] is True


def test_typer_expansion_outputs_json_and_latex(tmp_path):
    import sympy as s
    from typer.testing import CliRunner

    from hypra.cli import app
    from hypra.serialization import encode_expr

    eps = s.Symbol('eps')
    request = tmp_path / 'request.json'
    request.write_text(json.dumps({'expression': encode_expr(1/(1-eps)),
                                   'epsilon': encode_expr(eps), 'through': 2}))
    output, latex = tmp_path / 'result.json', tmp_path / 'result.tex'
    result = CliRunner().invoke(app, ['expand',str(request),'--output',str(output),'--latex',str(latex)])
    assert result.exit_code == 0, result.output
    assert json.loads(output.read_text())
    assert r'\[' in latex.read_text()


def test_typer_runtime_errors_still_json(tmp_path):
    from typer.testing import CliRunner

    from hypra.cli import app

    result = CliRunner().invoke(app, ['expand', str(tmp_path / 'missing.json')])
    assert result.exit_code == 2
    assert json.loads(result.stderr)['status'] == 'error'


def test_main_invalid_backend_preserves_usage_exit():
    import pytest

    with pytest.raises(SystemExit) as error:
        main(['backend-info','--backend','not-an-engine'])
    assert error.value.code == 2
