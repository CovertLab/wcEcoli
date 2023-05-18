"""Unit test for the chunk module."""

import io
from typing import cast, IO
import unittest

from wholecell.io.chunk import Chunk, FileLike


class Readonly:
    """Delegate to `stream` but support only `read()`; no `seek()` or `tell()`.
    """
    def __init__(self, stream: IO[bytes]) -> None:
        self.stream = stream

    def read(self, size: int = -1) -> bytes:
        return self.stream.read(size)


class Test_Chunk(unittest.TestCase):

    def test_read(self):
        """Test reading chunks, chunk type IDs, chunk sizes, and alignment."""
        f = io.BytesIO(
            b"MT01\x00\x00\x00\x00"
            b"DATA\x00\x00\x00\x03dat"
            b"0"
            b"MORE\x00\x00\x00\x04more")

        c1 = Chunk(f)
        assert c1.getname() == b'MT01'
        assert c1.getsize() == 0
        assert c1.read() == b''
        c1.close()

        c2 = Chunk(f)
        assert c2.getname() == b'DATA'
        assert c2.getsize() == 3
        b2 = c2.read()
        assert b2 == b'dat'  # the pad byte is not part of the data
        assert c2.read() == b''  # no more to read
        c2.close()

        c3 = Chunk(f)
        assert c3.getname() == b'MORE'
        assert c3.getsize() == 4
        b3 = c3.read()
        assert b3 == b'more'
        c3.close()

        with self.assertRaises(EOFError):
            Chunk(f)

        assert f.read() == b''  # at EOF

        if pow(0, 10):  # type-check this file-reading code but don't run it
            cf = Chunk(open('no_such_file.iff', 'rb'))
            cf.read()

    def test_container(self):
        """Test an IFF FORM containing a FORM-type code and 2 nested chunks.
        Also test seek() and tell().
        """
        CONTENT = (
            b"TEST"  # the FORM type
            b"MORE\x00\x00\x00\x04more"
            b"DATA\x00\x00\x00\x03dat"
            b"0")
        f = io.BytesIO(b"FORM\x00\x00\x00\x1c" + CONTENT)

        form = Chunk(f)
        assert form.getname() == b'FORM'
        assert form.getsize() == len(CONTENT)
        assert form.tell() == 0

        form_type = form.read(4)
        assert form_type == b'TEST'
        assert form.tell() == 4

        c1 = Chunk(form)
        assert c1.getname() == b'MORE'
        assert c1.read() == b'more'
        c1.close()
        assert form.tell() == 4 + (8 + 4)

        c2 = Chunk(form)
        assert c2.getname() == b'DATA'
        assert c2.read() == b'dat'
        assert c2.tell() == 3
        assert form.tell() == len(CONTENT)

        pos = c2.seek(0, 0)  # absolute seek c2 to the start of its chunk data
        assert pos == 0
        assert c2.tell() == 0
        assert form.tell() == len(CONTENT) - (3 + 1)

        assert c2.read() == b'dat'
        assert c2.tell() == 3
        assert form.tell() == len(CONTENT)

        pos = c2.seek(-3, 1)  # relative seek c2 to the start of its chunk data
        assert pos == 0
        assert c2.read() == b'dat'
        assert c2.tell() == 3
        assert form.tell() == len(CONTENT)  # including c2's pad byte
        c2.close()
        assert form.read() == b''  # at EOF

        pos = form.seek(-(8 + 3 + 1), 1)  # relative seek over c2 and its pad
        assert pos == 4 + (8 + 4)
        assert form.tell() == 4 + (8 + 4)

        c2a = Chunk(form)
        assert c2a.getname() == b'DATA'
        assert c2a.read() == b'dat'
        assert c2a.tell() == 3
        assert form.tell() == len(CONTENT)
        c2a.close()

        form.close()
        with self.assertRaises(ValueError):
            form.tell()

        assert f.read() == b''  # at EOF

    def test_options(self):
        """Test the Chunk options. For testing convenience, the variations are
        combined in a single test "file." Such mixed up contents make no sense
        in a real file since readers couldn't predict which options to use for
        each chunk. The options are unfortunate file format design variations.
        """
        f = io.BytesIO(
            b"UALN\x00\x00\x00\x03dat"  # unaligned, so no pad byte
            b"LITL\x07\x00\x00\x00little.0"  # little endian, aligned w/1 pad byte
            b"INCL\x0e\x00\x00\x00why???")  # includes the header in the chunk size

        c1 = Chunk(f, align=False)
        assert c1.getname() == b'UALN'
        assert c1.getsize() == 3
        b1 = c1.read()
        assert b1 == b'dat'
        c1.close()

        c2 = Chunk(f, bigendian=False)
        assert c2.getname() == b'LITL'
        assert c2.getsize() == 7
        b2 = c2.read()
        assert b2 == b'little.'
        c2.close()

        c3 = Chunk(f, bigendian=False, inclheader=True)
        assert c3.getname() == b'INCL'
        assert c3.getsize() == 6
        b3 = c3.read()
        assert b3 == b'why???'
        c3.close()

        with self.assertRaises(EOFError):
            Chunk(f)

        assert f.read() == b''  # at EOF

    def test_unseekable(self):
        """Test reading from a stream that doesn't support `seek()` or `tell()`.
        Chunks on such a stream can `tell()` but not `seek()`.
        """
        # Pretend this stream is FileLike for the type checker.
        # Chunk() will test if it supports `seek()` + `tell()`.
        f = cast(FileLike,
                 Readonly(io.BytesIO(
                    b"DATA\x00\x00\x00\x03dat"
                    b"0"
                    b"MORE\x00\x00\x00\x04more")))

        with self.assertRaises(AttributeError):
            f.tell()

        with self.assertRaises(AttributeError):
            f.seek(0, 0)

        c2 = Chunk(f)
        assert c2.getname() == b'DATA'
        assert c2.tell() == 0
        assert c2.read() == b'dat'
        assert c2.tell() == 3

        with self.assertRaises(OSError):
            c2.seek(0, 0)

        c2.close()
