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
