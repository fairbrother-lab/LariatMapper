COMP_NTS = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N'}
def reverse_complement(seq):
	return ''.join([COMP_NTS[seq[i]] for i in range(len(seq)-1,-1,-1)])


def comma_join(items) -> str:
	return ','.join(set(items))


def align_is_reverse(flag:int) -> bool:
	bit_flags = bin(int(flag))
	is_reverse = True if len(bit_flags)>=7 and bit_flags[-5]=='1' else False
	return is_reverse