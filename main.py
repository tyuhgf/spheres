import sys
from subprocess import Popen, PIPE
import threading
from time import sleep

p = Popen('/bin/bash', stdin=PIPE, stdout=PIPE, stderr=sys.stderr, shell=True)

stop = False


def write_all():
    while True:
        sleep(1)
        sys.stderr.write('123\n')
        data = p.stdout.read(1).decode("utf-8")
        if stop:
            print('STOP')
            break
        if data:
            sys.stderr.write(data)


writer = threading.Thread(target=write_all, args=())
writer.start()

while True:
    d = sys.stdin.readline().strip()
    if not d:
        print('terminate')
        stop = True
        sleep(2)
        writer.join(.5)
        break
    p.stdin.write((d + '\n').encode())
    p.stdin.flush()


#
#
#
# class LocalShell(object):
#     def __init__(self):
#         pass
#
#     def run(self):
#         env = os.environ.copy()
#         p = Popen('/bin/bash', stdin=PIPE, stdout=PIPE, stderr=subprocess.STDOUT, shell=True, env=env)
#         sys.stdout.write("Started Local Terminal...\r\n\r\n")
#
#         def writeall(p):
#             while True:
#                 # print("read data: ")
#                 data = p.stdout.read(1).decode("utf-8")
#                 if not data:
#                     break
#                 sys.stdout.write(data)
#                 sys.stdout.flush()
#
#         writer = threading.Thread(target=writeall, args=(p,))
#         writer.start()
#
#         try:
#             while True:
#                 d = sys.stdin.read(1)
#                 if not d:
#                     break
#                 self._write(p, d.encode())
#
#         except EOFError:
#             pass
#
#     def _write(self, process, message):
#         process.stdin.write(message)
#         process.stdin.flush()
#
#
# shell = LocalShell()
# shell.run()
