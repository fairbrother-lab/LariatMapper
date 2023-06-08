import os

shell = os.environ.get('SHELL')

if shell.endswith('bash'):
    bashrc_path = os.path.expanduser("~/.bashrc")
    with open(bashrc_path, "a") as bashrc:
        bashrc.write(
            "\nalias larmap='./larmap_top.sh'\n"
        )

elif shell.endswith('zsh'):
    zshrc_path = os.path.expanduser("~/.zshrc")
    with open(zshrc_path, "a") as zshrc:
        zshrc.write(
            "\nalias larmap='./larmap_top.sh'\n"
        )
