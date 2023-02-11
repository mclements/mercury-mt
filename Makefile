
test: test_mt.m mt.m
	mmc --make libmt test_mt && ./test_mt

PROXY += test
