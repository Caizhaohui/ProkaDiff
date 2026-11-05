# Parity Closure Supervisor

Independent orchestration for G1–G5 (locks, checkpoints, `agent -p` Controller rounds, Grok dispatch hooks).

## Verify CLI

```bash
python3 benchmark/results/parity_closure/supervisor/parity_closure_supervisor.py --verify-cli
```

## Tests (required before claiming auto-resume)

```bash
python3 benchmark/results/parity_closure/supervisor/parity_closure_supervisor.py --continuation-test
python3 benchmark/results/parity_closure/supervisor/parity_closure_supervisor.py --pause-test
```

## Run loop (login node: orchestration only; no cargo/WGS here)

```bash
nohup python3 benchmark/results/parity_closure/supervisor/parity_closure_supervisor.py \
  --loop --sleep 300 >> benchmark/results/parity_closure/state/supervisor_loop.log 2>&1 &
echo $! > benchmark/results/parity_closure/state/supervisor_loop.pid
```

## Operator control

```bash
benchmark/results/parity_closure/supervisor/parity_closure_ctl.sh pause|resume|stop|status
```

State: `benchmark/results/parity_closure/state/parity_closure_state.json`  
Control: `benchmark/results/parity_closure/state/supervisor_control.json`  
Floor anchor: `benchmark/results/parity_closure/state/baseline_w3_floors.json`
