from setuptools import find_packages, setup
import setuptools



setup(
        name='leafcutterITI',
        packages=find_packages(),
        version='0.1.8',
        description='LeafcutterITI implementation',
        author='Xingpei Zhang, David A Knowles',
        license='MIT',
        entry_points={
            'console_scripts': [
                'leafcutterITI-map = leafcutter.__main__:leafcutterITI_map_gen',
                'leafcutterITI-cluster = leafcutter.__main__:leafcutterITI_clustering'
            ],
        },
    )

