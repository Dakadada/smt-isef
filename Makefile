.PHONY: smoke-test test judge-demo
smoke-test test:
	PYTHONPATH=src /usr/bin/python3.12 -m unittest discover -s tests -v

judge-demo:
	test -n "$(RUN_DIR)" && test -n "$(PANEL)" && test -n "$(BUILD_TIMESTAMP)"
	PYTHONPATH=src /usr/bin/python3.12 scripts/judge_demo.py --run-dir "$(RUN_DIR)" --panel "$(PANEL)" --out "$${OUT:-artifacts/judge-demo/snapshot}" --build-timestamp "$(BUILD_TIMESTAMP)" $${ALLOW_INCOMPLETE:+--allow-incomplete-snapshot}
