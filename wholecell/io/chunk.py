"""Simple chunk reader for the EA IFF 85 file format and similar chunk formats.
Used in file formats such as AIFF, ILBM, RIFF, RMFF (RealMedia File Format),
and WAV.

© Copyright 2001-2023, Python Software Foundation.

Forked from the Python Standard Library's `chunk` module which is deprecated
since Python 3.11. This fork is a drop-in replacement. This v1.0 code is the
same except:
  * Chunk's `file` parameter is typed. It uses a `Protocol`, requiring
    Python 3.8+.
  * `seek()` returns the new position, per the `IO` protocol.
  * Bug fix: Relative seek now works after `read()` read the pad byte.
  * Clarified and corrected docstrings and comments to match the code and the
    file format standards.

See https://1fish2.github.io/IFF/
See https://github.com/1fish2/IFF
See http://fileformats.archiveteam.org/wiki/IFF

An IFF chunk has the following structure:

+----------------+
| ID (4 bytes)   |
+----------------+
| size (4 bytes) |
+----------------+
| data           |
| ...            |
+----------------+
[0]

The chunk type ID is a four character code that identifies the type of chunk,
represented as a 4-byte Python `bytes` object.

The chunk size is a 32-bit integer indicating the number of chunk data bytes.

If alignment is enabled, as in the IFF 85 standard, and the chunk size is odd,
then a 0 pad byte follows the chunk to align the next chunk to an even boundary.

An IFF-type file consists of one or more chunks.  In the IFF 85 standard, a file
consists of one container chunk.
"""

import struct
from typing import Protocol


class FileLike(Protocol):
    def read(self, size: int = -1) -> bytes: ...

    def seek(self, pos: int, whence: int = 0) -> int: ...

    def tell(self) -> int: ...


class Chunk:
    """
    A Chunk is a file-like object for reading an IFF chunk, implementing
        read(), close(), seek(), tell(), isatty().

    It also implements

        getname(): Return the chunk type ID -- a four character code as a bytes
            object.

        getsize(): Return the chunk size -- number of data bytes.

        skip(): Skip to the end of the chunk. Called by close().

    Use a `Chunk` instance per chunk within a file or container chunk.  At the
    end of the sequence, creating the next `Chunk` instance will raise EOFError.

    Example usage
        while True:
            try:
                chunk = Chunk(file)
            except EOFError:
                break

            chunktype = chunk.getname()
            data = chunk.read()
            # do something with data
            chunk.close()
    """
    closed: bool
    align: bool
    file: FileLike
    chunkname: bytes
    chunksize: int
    size_read: int
    offset: int
    seekable: bool

    def __init__(self, file: FileLike, align: bool = True,
                 bigendian: bool = True, inclheader: bool = False) -> None:
        """
        Instantate a Chunk to read an IFF chunk from a FileLike object such as
        a file or a container Chunk.

        Raises EOFError at the end of `file`.

        Args:
            file: The FileLike object to read from.  Actually, `read()` is the
              only required method.  If the `tell()` method is present and
              doesn’t raise an exception, Chunk will expect that `tell()`
              doesn't alter the object and that the `seek()` method is
              available.
            align: Whether chunks are aligned on 2-byte boundaries.
              Standard IFF 85 chunks are aligned.
            bigendian: Whether the chunk size field is stored in big endian
              format.  Standard IFF 85 chunks are big endian.  WAV and other
              RIFF file chunks are little endian.
            inclheader: Whether the chunk header (type and size) fields are
              counted in the chunk's size count, as in RMFF (RealMedia File
              Format) chunks; not in Standard IFF 85 chunks.
        """
        self.closed = False
        self.align = align
        if bigendian:
            strflag = '>'
        else:
            strflag = '<'
        self.file = file
        self.chunkname = file.read(4)
        if len(self.chunkname) < 4:
            raise EOFError
        try:
            self.chunksize = struct.unpack_from(strflag + 'L', file.read(4))[0]
        except struct.error:
            raise EOFError from None
        if inclheader:
            self.chunksize = self.chunksize - 8
        self.size_read = 0
        try:
            self.offset = self.file.tell()
        except (AttributeError, OSError):
            self.seekable = False
        else:
            self.seekable = True

    def getname(self) -> bytes:
        """Return the chunk type ID -- a four character code as a bytes object."""
        return self.chunkname

    def getsize(self) -> int:
        """Return the chunk size -- number of data bytes."""
        return self.chunksize

    def close(self) -> None:
        """Close and skip() to the end of the chunk.  This does not close the
        underlying file.
        """
        if not self.closed:
            try:
                self.skip()
            finally:
                self.closed = True

    def isatty(self) -> bool:
        """Returns False.  Raise ValueError if the close() method was called."""
        # NOTE: The Python Standard Library doc claimed that this and following
        #  methods raise OSError if the chunk is closed.  That would make more
        #  sense but it'd be an API change.
        if self.closed:
            raise ValueError("I/O operation on closed chunk")
        return False

    def seek(self, pos: int, whence: int = 0) -> int:
        """Seek to the specified position into the chunk.
        The initial position is 0 (the start of chunk).
        Raise OSError if the underlying file is not seekable.
        Raise ValueError if the close() method was called.

        Contrary to Python Standard Library doc, this method never implemented
        forward seeking over an unseekable file.

        Return the new absolute position.

        Args:
            pos: The target position into the chunk.
            whence: Defaults to 0 (absolute positioning); other values are 1
              (seek relative to the current position) and 2 (seek relative to
              the end of the chunk).
        """

        if self.closed:
            raise ValueError("I/O operation on a closed chunk")
        if not self.seekable:
            raise OSError("cannot seek")
        if whence == 1:
            pos = pos + min(self.size_read, self.chunksize)
        elif whence == 2:
            pos = pos + self.chunksize
        if pos < 0 or pos > self.chunksize:
            raise RuntimeError
        self.file.seek(self.offset + pos, 0)
        self.size_read = pos
        return pos

    def tell(self) -> int:
        """Return the current position into the chunk's data bytes.
        Raise ValueError if the close() method was called.
        """
        if self.closed:
            raise ValueError("I/O operation on a closed chunk")
        return min(self.size_read, self.chunksize)

    def read(self, size: int = -1) -> bytes:
        """Read at most `size` bytes from the chunk.
        If `size` is omitted or negative, read until the end of the chunk.
        Return an empty `bytes` object if already at the end of the chunk.
        Handle any pad byte for alignment.
        Raise ValueError if the close() method was called.
        """

        if self.closed:
            raise ValueError("I/O operation on a closed chunk")
        if self.size_read >= self.chunksize:
            return b''
        if size < 0:
            size = self.chunksize - self.size_read
        if size > self.chunksize - self.size_read:
            size = self.chunksize - self.size_read
        data = self.file.read(size)
        self.size_read = self.size_read + len(data)
        if (self.size_read == self.chunksize and
                self.align and
                (self.chunksize & 1)):
            dummy = self.file.read(1)
            self.size_read = self.size_read + len(dummy)
        return data

    def skip(self) -> None:
        """Skip the rest of the chunk, to position the underlying file to read
        the next chunk.  All further calls to read() will return `b''`.
        NOTE: close() will call skip() for you.  If you are not interested in
        the contents of the chunk, just call close().

        Raise ValueError if the close() method was called.
        """

        if self.closed:
            raise ValueError("I/O operation on a closed chunk")
        if self.seekable:
            try:
                n = self.chunksize - self.size_read
                if self.align and (self.chunksize & 1):
                    n = n + 1
                self.file.seek(n, 1)
                self.size_read = self.size_read + n
                return
            except OSError:
                pass
        while self.size_read < self.chunksize:
            n = min(8192, self.chunksize - self.size_read)
            dummy = self.read(n)
            if not dummy:
                raise EOFError
