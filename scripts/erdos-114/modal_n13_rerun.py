from pathlib import Path
import subprocess
import time

import modal


APP_NAME = "ehp114-n13-rerun-20260507"
VOLUME_NAME = "ehp114-n13-rerun-20260507"
REMOTE_SRC = "/workspace/erdos-114"
REMOTE_OUT = "/outputs/EXP-MM-EHP-007-n13-inari-MODAL-20260507-01"

ROOT = Path(__file__).resolve().parent

volume = modal.Volume.from_name(VOLUME_NAME, create_if_missing=True)

image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install(
        "ca-certificates",
        "curl",
        "build-essential",
        "pkg-config",
        "m4",
        "autoconf",
        "automake",
        "libtool",
        "make",
    )
    .run_commands(
        "curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | "
        "sh -s -- -y --profile minimal"
    )
    .env(
        {
            "PATH": "/root/.cargo/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin",
            "RUSTFLAGS": "-Ctarget-cpu=haswell",
        }
    )
    .add_local_file(ROOT / "Cargo.toml", f"{REMOTE_SRC}/Cargo.toml", copy=True)
    .add_local_file(ROOT / "Cargo.lock", f"{REMOTE_SRC}/Cargo.lock", copy=True)
    .add_local_dir(ROOT / "src", f"{REMOTE_SRC}/src", copy=True)
    .run_commands(f"cd {REMOTE_SRC} && cargo build --release --bin ehp_general_ieee1788")
)

app = modal.App(APP_NAME)


def _run_streaming(cmd: list[str], commit_output_volume: bool = False) -> dict:
    print("MODAL_CMD", " ".join(cmd), flush=True)
    start = time.time()
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )
    assert proc.stdout is not None
    for line in proc.stdout:
        print(line.rstrip(), flush=True)
        if commit_output_volume:
            try:
                volume.commit()
            except Exception as exc:
                print(f"WARNING volume commit failed during stream: {exc}", flush=True)
    rc = proc.wait()
    if commit_output_volume:
        try:
            volume.commit()
        except Exception as exc:
            print(f"WARNING final volume commit failed: {exc}", flush=True)
    return {"returncode": rc, "elapsed_secs": time.time() - start}


@app.function(image=image, cpu=4, memory=8192, timeout=30 * 60)
def smoke_modal() -> dict:
    outdir = Path("/tmp/ehp114-modal-smoke-n3")
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [
        f"{REMOTE_SRC}/target/release/ehp_general_ieee1788",
        "--outdir",
        str(outdir),
        "--only",
        "3",
        "--checkpoint-every",
        "512",
        "--log-every",
        "512",
    ]
    result = _run_streaming(cmd)
    result["files"] = sorted(str(p.relative_to(outdir)) for p in outdir.glob("*"))
    if result["returncode"] != 0:
        raise RuntimeError(f"smoke failed with return code {result['returncode']}")
    return result


@app.function(
    image=image,
    volumes={"/outputs": volume},
    cpu=32,
    memory=65536,
    timeout=24 * 60 * 60,
    nonpreemptible=True,
)
def run_n13_modal(
    resume: bool = True,
    checkpoint_every: int = 25_000,
    log_every: int = 25_000,
) -> dict:
    outdir = Path(REMOTE_OUT)
    outdir.mkdir(parents=True, exist_ok=True)

    cmd = [
        f"{REMOTE_SRC}/target/release/ehp_general_ieee1788",
        "--outdir",
        str(outdir),
        "--only",
        "13",
        "--checkpoint-every",
        str(checkpoint_every),
        "--log-every",
        str(log_every),
    ]
    if resume:
        cmd.append("--resume")

    print("MODAL_N13_CMD", " ".join(cmd), flush=True)
    result = _run_streaming(cmd, commit_output_volume=True)
    files = sorted(str(p.relative_to(outdir)) for p in outdir.glob("*"))
    result["outdir"] = str(outdir)
    result["files"] = files
    return result


@app.local_entrypoint()
def smoke(wait: bool = True):
    call = smoke_modal.spawn()
    print(
        {
            "call_id": call.object_id,
            "dashboard_url": call.get_dashboard_url(),
            "wait": wait,
        }
    )
    if wait:
        print(call.get(timeout=30 * 60))


@app.local_entrypoint()
def main(
    resume: bool = True,
    checkpoint_every: int = 25_000,
    log_every: int = 25_000,
    wait: bool = False,
):
    call = run_n13_modal.spawn(
        resume=resume,
        checkpoint_every=checkpoint_every,
        log_every=log_every,
    )
    print(
        {
            "call_id": call.object_id,
            "dashboard_url": call.get_dashboard_url(),
            "resume": resume,
            "checkpoint_every": checkpoint_every,
            "log_every": log_every,
            "wait": wait,
        }
    )
    if wait:
        print(call.get(timeout=24 * 60 * 60))
