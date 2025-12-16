import yaml
from pathlib import Path
from qss.distributions import DistType
from qss.solver_options import CoolingSchedule

def load_config(config_name):
    """Loads a YAML config file from the configs directory."""
    config_path = Path(__file__).parent.parent / "config" / f"{config_name}.yaml"
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    # Process Enums
    if "physics" in config and "dist_type" in config["physics"]:
        dist_str = config["physics"]["dist_type"]
        config["physics"]["dist_type"] = getattr(DistType, dist_str)

    if "solver" in config:
        # Handle CoolingSchedule
        if "sched" in config["solver"]:
            sched_str = config["solver"]["sched"]
            config["solver"]["sched"] = getattr(CoolingSchedule, sched_str)

        # Handle bounds (list of lists -> list of tuples)
        if "bounds" in config["solver"]:
            config["solver"]["bounds"] = [tuple(b) for b in config["solver"]["bounds"]]

    return config




if __name__ == "__main__":
    sample_config = load_config("")
    print(sample_config)
