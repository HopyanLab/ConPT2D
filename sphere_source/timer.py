#!/usr/bin/env python3
import time

################################################################################

class timer_error(Exception):
	"""custom exception used to report errors in use of timer class"""

################################################################################

class timer:
	def __init__(self):
		self._start_time = None
	def start(self):
		"""start a new timer"""
		if self._start_time is not None:
			raise timer_error(f"timer is running, use .stop()")
		self._start_time = time.perf_counter()
	def stop(self):
		"""stop the timer, and report the elapsed time"""
		if self._start_time is None:
			raise timer_error(f"timer is not running, use .start()")
		elapsed_time = time.perf_counter() - self._start_time
		self._start_time = None
		hours = int(elapsed_time/3600)
		minutes = int((elapsed_time - hours*3600)/60)
		seconds = elapsed_time - hours*3600 - minutes*60
		print(f"elapsed time: {hours:d}:{minutes:02d}:{seconds:02.4f}")

################################################################################

if __name__ == '__main__':
	pass

################################################################################
# EOF
