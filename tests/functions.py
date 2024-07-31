import gzip

def files_equal(file_a, file_b, sort=False):
	# Load
	if str(file_a).endswith('.gz'):
		with gzip.open(file_a, 'rt') as in_file:
			lines_a = in_file.readlines()
	else:
		with open(file_a) as in_file:
			lines_a = in_file.readlines()
			
	if str(file_b).endswith('.gz'):
		with gzip.open(file_b, 'rt') as in_file:
			lines_b = in_file.readlines()
	else:
		with open(file_b) as in_file:
			lines_b = in_file.readlines()

	# Count
	linecount = len(lines_a)
	assert linecount == len(lines_b)

	# Sort
	if sort is True:
		lines_a.sort()
		lines_b.sort()

	# Compare
	assert lines_a == lines_b, f'{file_a} vs {file_b}'