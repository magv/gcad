from setuptools import setup
from setuptools.command.build_ext import build_ext
import subprocess

class my_build_ext(build_ext):
    def run(self):
        subprocess.check_call(["make", "build-deps"])
        super().run()

setup(cmdclass={"build_ext": my_build_ext})
