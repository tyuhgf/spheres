import asyncio
from settings import gap_path


class Communicator:
    def __init__(self, process):
        self.counter = 0
        self.process = process

        self.answers = ''
        self.questions = []
        self._step = .3

        self.listener = None

    async def process_errors(self):
        while True:
            try:
                if self.process.returncode is not None:
                    break
                line = await asyncio.wait_for(self.process.stderr.read(200), 1)
                line = line.decode(errors='ignore')
                if line:
                    print('GAP ERROR > ', line)
                    self.listener.error(line)
            except asyncio.TimeoutError:
                pass
        pass

    async def process_answers(self):
        while True:
            try:
                if self.process.returncode is not None:
                    break
                line = await asyncio.wait_for(self.process.stdout.read(200), 1)
                line = line.decode(errors='ignore')
                self.answers += line
                if self.listener is not None:
                    self.listener(line)
            except asyncio.TimeoutError:
                pass
        pass

    async def request(self, lines, timeout=1, listener=None):
        if not isinstance(lines, list):
            lines = [lines]
        self.listener = listener
        n_answers = len(self.answers)
        for line in lines:
            if self.process.returncode is not None:
                return None, 'terminated', {'returncode': self.process.returncode}
            self.process.stdin.write((line + '\n').encode())
        for _ in range(int(timeout // self._step + 1)):
            await asyncio.sleep(self._step)
            if self.listener is not None and self.listener.complete:
                res = self.answers[n_answers:].split('\x1b[1m\x1b[34mgap> \x1b[0m\x1b[31m')  # gap prompt string
                if res[-1] == '':
                    res = res[:-1]
                return res, 'ok', {}
        return self.answers[n_answers:], 'ok', {'timeout': True}


class Listener:
    def __init__(self):
        self.complete = False

    def __call__(self, _line):
        self.complete = True

    def error(self, _line):
        self.complete = True


class ListenerNLines(Listener):
    def __init__(self, n=1):
        super(ListenerNLines, self).__init__()
        self.n = n
        if n == 0:
            self.complete = True
        self.answers = ''

    def __call__(self, line):
        self.answers += line
        if self.answers.count('gap> ') >= self.n:
            self.complete = True


async def execute_series(requests: list):
    if not isinstance(gap_path, str) and not isinstance(gap_path, bytes):
        raise ValueError('Invalid path to GAP executable!')

    loop = asyncio.get_event_loop()

    p = await asyncio.create_subprocess_exec(
        gap_path,
        stdin=asyncio.subprocess.PIPE,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE)

    c = Communicator(p)

    loop.create_task(c.process_answers())
    loop.create_task(c.process_errors())
    await asyncio.sleep(5)  # to start gap system

    answers = []
    for req in requests:
        if isinstance(req, str):
            ans = await loop.create_task(c.request(req, listener=ListenerNLines()))
        elif isinstance(req, tuple):
            ans = await loop.create_task(c.request(*req, listener=ListenerNLines()))
        elif isinstance(req, dict):
            ans = await loop.create_task(c.request(**req))
        else:
            raise ValueError('Unknown type of request!')

        if ans[1] != 'ok':
            raise Exception(f'Gap communication error {ans}')
        answers.append(ans)

    p.terminate()
    return answers


def gap_execute_commands(requests):
    return asyncio.get_event_loop().run_until_complete(execute_series(requests))
