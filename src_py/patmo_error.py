import sys
import inspect


# error trigger, requires fileName as __file__
# and frame as inspect.currentframe()
# Example call is
#   patmo_error.trigError(__file__,inspect.currentframe())
def trigError(file_name, frame):
    src_file_name = file_name.replace(".pyc", ".py")
    print "(Triggered by source file: " + src_file_name + " [" \
          + str(inspect.getlineno(frame)) + "])"
    # stop program execution
    sys.exit()
