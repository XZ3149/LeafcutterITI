from setuptools import find_packages, setup
import setuptools



setup(
        name='tealeaf',
        packages=find_packages(),
        version='0.1.0',
        description='tealeaf implementation',
        author='Xingpei Zhang, David A Knowles',
        license='MIT',
        entry_points={
            'console_scripts': [
                'tealeaf-map = tealeaf.__main__:tealeaf_map_gen',
                'tealeaf-cluster = tealeaf.__main__:tealeaf_clustering',
                'tealeaf-sc = tealeaf.__main__:tealeaf_sc',
                'tealeaf-ggsashimi = tealeaf.__main__:tealeaf_ggsashimi'
            ],
        },
    )

