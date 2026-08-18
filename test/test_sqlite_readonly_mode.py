from singlem.metapackage_read_name_store import MetapackageReadNameStore
from singlem.sequence_database import SequenceDatabase


def test_sequence_database_acquire_uses_readonly_immutable_sqlite(monkeypatch, tmp_path):
    db_dir = tmp_path / "db"
    db_dir.mkdir()
    (db_dir / "CONTENTS.json").write_text('{"singlem_database_version": 5}')
    (db_dir / "otus.sqlite3").touch()

    captured = {}

    class FakeEngine:
        def connect(self):
            return object()

    def fake_create_engine(url):
        captured["url"] = url
        return FakeEngine()

    monkeypatch.setattr("singlem.sequence_database.create_engine", fake_create_engine)

    SequenceDatabase.acquire(str(db_dir))

    assert captured["url"].startswith("sqlite+pysqlite:///file:")
    assert "mode=ro" in captured["url"]
    assert "immutable=1" in captured["url"]
    assert "uri=true" in captured["url"]


def test_metapackage_read_name_store_acquire_uses_readonly_immutable_sqlite(monkeypatch, tmp_path):
    sqlite_db = tmp_path / "read_taxonomies.sqlite3"
    sqlite_db.touch()

    captured = {}

    def fake_create_engine(url, echo=None, future=None):
        captured["url"] = url
        captured["echo"] = echo
        captured["future"] = future
        return object()

    monkeypatch.setattr("singlem.metapackage_read_name_store.create_engine", fake_create_engine)

    MetapackageReadNameStore.acquire(str(sqlite_db))

    assert captured["url"].startswith("sqlite+pysqlite:///file:")
    assert "mode=ro" in captured["url"]
    assert "immutable=1" in captured["url"]
    assert "uri=true" in captured["url"]
    assert captured["future"] is True
