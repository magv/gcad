from setuptools import setup
from setuptools.command.build_ext import build_ext
import subprocess
import sys

class my_build_ext(build_ext):
    def run(self):
        subprocess.check_call(["make", "build-deps"])
        super().run()

    def build_extensions(self):
        if sys.platform == "darwin":
            extra_link_args = ["-Wl,-dead_strip"]
        else:
            extra_link_args = ["-Wl,--gc-sections"]
        for ext in self.extensions:
            ext.extra_link_args.extend(extra_link_args)
        super().build_extensions()

setup(cmdclass={"build_ext": my_build_ext})
