from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_nginx_exposes_biostar_as_a_controller_prefix() -> None:
    configuration = (
        PROJECT_ROOT / "nginx" / "default.conf.template"
    ).read_text(encoding="utf-8")

    assert "server_name ${BIOSTAR_DOMAIN};" in configuration
    assert "location = /biostar {" in configuration
    assert "location /biostar/ {" in configuration
    assert "proxy_pass http://api:8000/;" in configuration
    assert "proxy_set_header X-Forwarded-Prefix /biostar;" in configuration
    assert "location / {" not in configuration
    assert "location /api/ {" not in configuration


def test_production_api_uses_biostar_as_root_path() -> None:
    compose = (
        PROJECT_ROOT / "docker-compose.prod.yml"
    ).read_text(encoding="utf-8")

    assert '"--root-path", "/biostar"' in compose


def test_production_stack_uses_only_sqlite() -> None:
    compose = (
        PROJECT_ROOT / "docker-compose.prod.yml"
    ).read_text(encoding="utf-8")

    assert "sqlite:////data/biostar.db" in compose
    assert "postgres" not in compose.lower()
    assert "mysql" not in compose.lower()
    assert "pgadmin" not in compose.lower()


def test_development_stack_uses_only_sqlite() -> None:
    compose = (
        PROJECT_ROOT / "docker-compose.dev.yml"
    ).read_text(encoding="utf-8")

    assert "sqlite:////data/biostar.db" in compose
    assert "postgres" not in compose.lower()
    assert "mysql" not in compose.lower()
    assert "pgadmin" not in compose.lower()


def test_environment_examples_do_not_define_external_database_services() -> None:
    for filename in (".env.dev.example", ".env.prod.example"):
        content = (PROJECT_ROOT / filename).read_text(encoding="utf-8").lower()

        assert "postgres" not in content
        assert "mysql" not in content
        assert "pgadmin" not in content


def test_production_stack_does_not_depend_on_removed_seed_service() -> None:
    compose = (
        PROJECT_ROOT / "docker-compose.prod.yml"
    ).read_text(encoding="utf-8")

    assert "\n  seed:" not in compose
    assert "condition: service_completed_successfully" in compose
    assert "depends_on:\n      migrate:" in compose


def test_development_stack_does_not_depend_on_removed_seed_service() -> None:
    compose = (
        PROJECT_ROOT / "docker-compose.dev.yml"
    ).read_text(encoding="utf-8")

    assert "\n  seed:" not in compose
    assert "depends_on:\n      migrate:" in compose
