# Parse Analysis Methodology Index

Each active or legacy workflow script in `analysis/` should have a matching methodology file under the corresponding subfolder here.

Required methodology content:

- Script purpose and biological/computational claim
- Exact inputs, external references, and downloads
- Processing steps in implemented order
- Cache/intermediate reuse behavior
- Output tiers: `intermediate/`, `tables/`, `figures/`, `logs/`, `reports/`
- Downstream dependencies and whether the script is terminal
- Current status: active, legacy comparison, or delete candidate

Shared project conventions are documented in:

- `analysis/script_inventory_and_dependency_map.md`
- `analysis/methodology/common/shared_configuration_and_logging_methodology.md`
