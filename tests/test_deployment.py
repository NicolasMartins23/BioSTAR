from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_nginx_exposes_biostar_as_a_controller_prefix() -> None:
    configuration = (
        PROJECT_ROOT / "nginx" / "default.conf.template"
    ).read_text(encoding="utf-8")

    assert "server_name ${BIOSTAR_DOMAIN};" in configuration
    assert "location = /biostar {" in configuration
    assert "location /biostar/ {" in configuration
    assert "proxy_pass http://api:8000;" in configuration
    assert "location / {" not in configuration
    assert "location /api/ {" not in configuration


def test_production_stack_does_not_depend_on_removed_seed_service() -> None:
    compose = (
        PROJECT_ROOT / "docker-compose.prod.yml"
    ).read_text(encoding="utf-8")

    assert "
  seed:" not in compose
    assert "condition: service_completed_successfully" in compose
    assert "depends_on:
      migrate:" in compose


def test_development_stack_does_not_depend_on_removed_seed_service() -> None:
    compose = (
        PROJECT_ROOT / "docker-compose.dev.yml"
    ).read_text(encoding="utf-8")

    assert "
  seed:" not in compose
    assert "depends_on:
      migrate:" in compose
