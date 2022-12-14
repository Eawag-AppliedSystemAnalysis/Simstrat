import csv
import math
import os
import shutil
import subprocess
import sys

INPUT_DATA_FOLDER = "TestCase_LakeZurich"
ORIGINAL_CONFIG_FILE = "TestCase_LakeZurich.par"
USED_CONFIG_FILE = "config.par"
SIMSTRAT_EXE = "../build/simstrat"
EXPECTED_RESULTS_FOLDER = "TestCases_Results_expected_windows" if sys.platform.startswith('win') else "TestCases_Results_expected"
ACTUAL_RESULTS_FOLDER = "TestCases_Results"
ERROR_FOLDER = "TestDiffs"

def read_lines(file):
    with open(file) as f:
        return f.readlines()

def write_lines(file, lines):
    with open(file, 'w') as f:
        for line in lines:
            f.write(line)

def as_float_or_null(str):
    try:
        return float(str)
    except Exception:
        return 0

def as_row(line):
    return [as_float_or_null(cell) for cell in line.split(",")]

def is_equal(x, y, tol=1e-8):
    if math.isnan(x):
        return math.isnan(y)
    else:
        min, max = (x, y) if abs(x) < abs(y) else (y,x)
        return abs(min/max - 1) < tol

def join_paths(prefix, path):
    return "%s/%s" % (prefix, path)

class Data(object):
    def __init__(self, name):
        self.name = name
        self.header = None
        self.rows = []

    def handle_row(self, row):
        if self.header is None:
            self.header = row
        else:
            self.rows.append([float(c) for c in row])

def read_data(file_path):
    data = Data(file_path)
    with open(file_path, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            data.handle_row(row)
    return data

class Tester(object):
    def __init__(self):
        shutil.rmtree(ERROR_FOLDER, ignore_errors=True)
        self.number_of_failures = 0

    def close(self):
        if self.number_of_failures > 0:
            print("Test FAILED: %s failures occurred" % self.number_of_failures)
            return 1
        print("Test SUCCESS")
        return 0

    def input_file(self, file):
        return join_paths(INPUT_DATA_FOLDER, file)
    
    def output_file(self, variable):
        return join_paths(ACTUAL_RESULTS_FOLDER, "%s_out.dat" % variable)

    def log_title(self, title):
        print()
        print("::::::::::::::::::::::::::::::::::::::::::")
        print(title)
        print("::::::::::::::::::::::::::::::::::::::::::")

    def log_step(self, msg):
        print("======= %s =======" % msg)

    def fail(self, msg):
        print("ERROR: %s" % msg)
        self.number_of_failures += 1

    def clean_up(self):
        self.log_step("Clean up")
        for f in os.listdir(ACTUAL_RESULTS_FOLDER):
            if f.endswith(".dat"):
                os.remove(join_paths(ACTUAL_RESULTS_FOLDER, f))
        if os.path.exists(self.input_file("forcing.dat")):
            os.remove(self.input_file("forcing.dat"))

    def extract_forcing_data_subset(self, end_time):
        lines = []
        with open(self.input_file("MeteoForcing.dat")) as f:
            for line in f:
                if len(lines) == 0 or float(line.split()[0]) <= end_time:
                    lines.append(line)
        write_lines(self.input_file("forcing.dat"), lines)

    def create_simulation_config(self, parameters):
        lines = []
        section = None
        with open(ORIGINAL_CONFIG_FILE) as f:
            for line in f:
                splitted_line = line.split(":", maxsplit=1)
                if len(splitted_line) == 1:
                    lines.append(line)
                else:
                    parameter = splitted_line[0].strip()[1:-1]
                    value = splitted_line[1].split("!")[0].strip()
                    if value == '{':
                        section = parameter
                    key = "%s.%s" % (section, parameter)
                    if key in parameters:
                        new_value = parameters[key]
                        if isinstance(new_value, str):
                            new_value = '"%s"' % new_value
                        elif new_value == True:
                            new_value = "true"
                        elif new_value == False:
                            new_value = "false"
                        eol = "," if value.endswith(",") else ""
                        lines.append('"%s":%s%s\n' % (parameter, new_value, eol))
                    else:
                        lines.append(line)
        write_lines(USED_CONFIG_FILE, lines)

    def run_simulation(self, end_time):
        self.log_step("Run until %s" % end_time)
        self.extract_forcing_data_subset(end_time)
        self.create_simulation_config({"Input.Forcing" : self.input_file("forcing.dat"), 
                                       "Output.All" : True,
                                       "Simulation.End d" : end_time})
        self.execute_simstrat(USED_CONFIG_FILE)

    def run_simulation_with_defined_times(self, times):
        self.log_step("Run with defined times %s" % times)
        self.create_simulation_config({"Output.Times" : times, "Output.All" : True})
        self.execute_simstrat(USED_CONFIG_FILE)

    def execute_simstrat(self, config_file):
        subprocess.run([SIMSTRAT_EXE, config_file], check=True)

    def read_line(self, file, time):
        lines = read_lines(file)
        for i in range(1, len(lines)):
            line = lines[i]
            if is_equal(as_float_or_null(line.split(",")[0]), time):
                return line
        raise ValueError("%s not found in %s" % (time, file))

    def read_row(self, file, time):
        return as_row(self.read_line(file, time))

    def read_and_interpolate(self, variable):
        path = self.output_file(variable)
        lines = read_lines(path)
        if len(lines) != 3:
            self.fail("Output file %s has %s lines instead of 3" % (path, len(lines)))
        r1 = as_row(lines[1])
        r0 = as_row(lines[2])
        return [0.2 * r1[i] + 0.8 * r0[i] for i in range(len(r0))]

    def assert_equal_outputs(self):
        for f in os.listdir(EXPECTED_RESULTS_FOLDER):
            if f.endswith("_out.dat"):
                actual_file = join_paths(ACTUAL_RESULTS_FOLDER, f)
                if os.path.exists(actual_file):
                    self.assert_almost_equal_result(f)
                else:
                    self.fail("Output file %s missing" % actual_file)

    def assert_almost_equal_result(self, file):
        expected_result_file = join_paths(EXPECTED_RESULTS_FOLDER, file)
        expected_result = read_data(expected_result_file)
        actual_result_file = join_paths(ACTUAL_RESULTS_FOLDER, file)
        actual_result = read_data(actual_result_file)
        if actual_result.header != expected_result.header:
            self.fail("Actual header in %s differs from expected header in %s."
                      % (actual_result_file, expected_result_file))
        elif len(actual_result.rows) != len(expected_result.rows):
            self.fail("Actual file %s has %s rows but expected file %s has % rows."
                      % (actual_result_file, len(actual_result.rows), expected_result_file, len(expected_result.rows)))
        else:
            for i in range(len(actual_result.rows)):
                actual_row = actual_result.rows[i]
                expected_row = expected_result.rows[i]
                if len(actual_row) != len(expected_row):
                    self.fail("Actual row at %s in %s has %s cells but expected row at %s in %s has %s cells."
                              % (actual_row[0], actual_result_file, len(actual_row), 
                                 expected_row[0], expected_result_file, len(expected_row)))
                    break
                if actual_row[0] != expected_row[0]:
                    self.fail("Actual row at %s in %s but expected row at %s in %s."
                              % (actual_row[0], actual_result_file, expected_row[0], expected_result_file))
                    break
                max_absolute_diff = 0
                max_relative_diff = 0
                for j in range(1, len(actual_row)):
                    ac = actual_row[j]
                    ec = expected_row[j]
                    absolute_diff = abs(ac - ec)
                    max_absolute_diff = max(max_absolute_diff, absolute_diff)
                    if absolute_diff > 0:
                        max_relative_diff = max(max_relative_diff, absolute_diff / max(abs(ac), abs(ec)))
                if max_relative_diff > 1e-3:
                    print("%s @ %.3f: %.2e, %.2e" % (file, actual_row[0], max_absolute_diff, max_relative_diff))
                    if max_absolute_diff > 1e-10:
                        self.fail("%s @ %.3f: max error = %.2e, max rel error = %.2e" 
                                  % (file, actual_row[0], max_absolute_diff, max_relative_diff))

    def assert_equal_lists(self, actual_list, expected_list, list_name):
        if len(actual_list) != len(expected_list):
            self.fail("%s: expected size %s != actual size %s" % (list_name, len(expect_list), len(actual_list)))
        else:
            for i in range(len(actual_list)):
                if not is_equal(actual_list[i], expected_list[i], 1e-4):
                    self.fail("%s: expected\n%s\nbut was\n%s" % (list_name, expected_list, actual_list))
                    break

    def assert_equal_outputs_for(self, time):
        self.log_step("Check output for time %s" % time)
        for f in os.listdir(EXPECTED_RESULTS_FOLDER):
            if f.endswith("_out.dat"):
                actual_file = join_paths(ACTUAL_RESULTS_FOLDER, f)
                if os.path.exists(actual_file):
                    expected_line = self.read_line(join_paths(EXPECTED_RESULTS_FOLDER, f), time)
                    actual_line = self.read_line(join_paths(ACTUAL_RESULTS_FOLDER, f), time)
                    if actual_line != expected_line:
                        self.fail("%s: expected\n%s\nbut was\n%s" % (f, expected_line, actual_line))
                else:
                    self.fail("Output file %s missing" % actual_file)


tester = Tester()

tester.log_title("Test incremental simulation with output at equidistant time steps")
tester.clean_up()
tester.run_simulation(12150.0)
tester.run_simulation(12250.125)
tester.run_simulation(12350.375)
tester.log_step("Run until end (12732)")
tester.create_simulation_config({"Output.All" : True})
tester.execute_simstrat(USED_CONFIG_FILE)
tester.assert_equal_outputs()

tester.log_title("Test incremental simulation and output interpolation for a list of output time points")
tester.clean_up()
tester.run_simulation_with_defined_times([12010.09722222222222, 12010.1006944444444])
v_expected = tester.read_and_interpolate("V")
nn_expected = tester.read_and_interpolate("NN")
hk_expected = tester.read_and_interpolate("HK")
tester.clean_up()
tester.run_simulation_with_defined_times([12002.5, 12002.625, 12003, 12010.09722222])
tester.run_simulation_with_defined_times([12002.5, 12002.625, 12003, 12010.09722222, 12010.1, 12010.1006944, 12020, 12100])
tester.log_step("Check interpolation")
tester.assert_equal_lists(tester.read_row(tester.output_file("V"), 12010.1), v_expected, "V")
tester.assert_equal_lists(tester.read_row(tester.output_file("NN"), 12010.1), nn_expected, "NN")
tester.assert_equal_lists(tester.read_row(tester.output_file("HK"), 12010.1), hk_expected, "HK")
tester.assert_equal_outputs_for(12003)
tester.assert_equal_outputs_for(12020)
tester.assert_equal_outputs_for(12100)

sys.exit(tester.close())