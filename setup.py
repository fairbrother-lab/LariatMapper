import yaml
from setuptools import setup, find_packages

def get_dependencies_from_yaml(yaml_file):
    with open(yaml_file, 'r') as f:
        env_data = yaml.safe_load(f)
    dependencies = env_data.get('dependencies', [])
    return dependencies


setup(
    name='larmap-test',
    version='1.0.0',
    packages=find_packages(),
    install_requires=get_dependencies_from_yaml('env.yaml'),
    author='nealneal',
    description='lariat mapping (test package)',
    url='https://github.com/nyin01/lariat_mapping_test',
)