from setuptools import find_packages, setup
import setuptools



setup(
        name='leafcutterITI',
        packages=find_packages(),
        version='0.2.0',
        description='LeafcutterITI implementation',
        author='Xingpei Zhang, David A Knowles',
        license='MIT',
        entry_points={
            'console_scripts': [
                'leafcutterITI-map = leafcutterITI.__main__:leafcutterITI_map_gen',
                'leafcutterITI-cluster = leafcutterITI.__main__:leafcutterITI_clustering',
                'leafcutterITI-scITI = leafcutterITI.__main__:leafcutterITI_scITI'
            ],
        },
    )

