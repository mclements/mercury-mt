
test: test_mt.m mt.m
	mmc --make test_mt && ./test_mt

clean:
	rm -rf Mercury
	rm -f test_mt_1 test_mt_2

PROXY += test
PROXY += clean
