# Run a script and report its peak resident memory -- /usr/bin/time is not installed here.
import resource, runpy, sys, time
t = time.time()
script = sys.argv[1]; sys.argv = sys.argv[1:]
try:
    runpy.run_path(script, run_name="__main__")
finally:
    print(f"[memrun] peak RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024:.0f} MB, "
          f"wall {(time.time()-t)/60:.1f} min", flush=True)
