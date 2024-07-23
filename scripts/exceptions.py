import subprocess

class RunCommandError(Exception):
	'''
	
	'''

	def __init__(self, process:subprocess.CompletedProcess):
		super().__init__()
		self.process = process
		self.response = process.stdout + process.stderr

	def __str__(self):
		return f'Command returned non-zero exit status {self.process.returncode}. \n{self.response}'
