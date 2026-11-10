from core.paths import PROJECT
import yaml

WORK = PROJECT


def load_config() -> dict:
    return yaml.safe_load((PROJECT / 'config/agency_validation.yaml').read_text(encoding='utf-8'))
