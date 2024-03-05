from setuptools import find_packages, setup

import setuptools

if __name__ == "__main__":
    setuptools.setup(py_modules = ["leafcutter"])

if False: 
    setup(
        name='leafcutterITI',
        packages=find_packages(),
        version='0.1.0',
        description='LeafcutterITI implementation',
        author='Xingpei Zhang, David A Knowles',
        license='MIT',
        entry_points={
            'console_scripts': [
                'leafcutter-ds = leafcutter.__main__:leafcutter_ds',
            ],
        },
    )

