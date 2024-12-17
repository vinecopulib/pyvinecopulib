"""
Serves documentation for `pyvinecopulib`.
"""

import argparse
import os
import sys
import webbrowser
from http.server import SimpleHTTPRequestHandler
from os.path import abspath, dirname, isdir, join
from shutil import rmtree
from socketserver import TCPServer
from subprocess import check_call

from gen_sphinx import write_doc_modules, write_examples


class _Handler(SimpleHTTPRequestHandler):
  # An HTTP handler without logging.
  def log_request(*_):
    pass


def _str2bool(value):
  # From: https://stackoverflow.com/a/19233287/7829525

  return value.lower() in ("yes", "y", "true", "t", "1")


def preview_main(gen_script, default_port):
  """Main entrypoint for previewing documentation.
  Args:
      gen_script: Generation script, required to generate docs.
      default_port: Default port for local HTTP server.
  """
  parser = argparse.ArgumentParser()
  parser.register("type", "bool", _str2bool)
  parser.add_argument(
    "--browser",
    type="bool",
    default=True,
    metavar="BOOL",
    help="Open browser. Disable this if you are frequently recompiling.",
  )
  parser.add_argument(
    "--port",
    type=int,
    default=default_port,
    metavar="PORT",
    help="Port for serving doc pages with a HTTP server.",
  )
  parser.add_argument(
    "--generates",
    type="bool",
    default=False,
    metavar="BOOL",
    help="Only generates.",
  )
  args = parser.parse_args()

  if args.generates:
    write_doc_modules(dirname(gen_script))
    write_examples(dirname(gen_script))
  else:
    # Choose an arbitrary location for generating documentation.
    out_dir = abspath("_build")

    if isdir(out_dir):
      rmtree(out_dir)
    # Generate.
    check_call([sys.executable, gen_script, "--out_dir", out_dir])
    print("Sphinx preview docs are available at:")
    file_url = "file://{}".format(join(out_dir, "index.html"))

    browser_url = file_url
    print()
    print("  {}".format(file_url))
    # Serve the current directory for local browsing. Required for MacOS.
    # N.B. We serve the preview via a HTTP server because it is necessary for
    # certain browsers (Safari on MacOS, possibly Chrome) due to local file
    # restrictions.
    os.chdir(out_dir)
    sockaddr = ("127.0.0.1", args.port)
    TCPServer.allow_reuse_address = True
    httpd = TCPServer(sockaddr, _Handler)
    http_url = "http://{}:{}/index.html".format(*sockaddr)
    print()
    print("  {}".format(http_url))
    # Default to using HTTP serving only on MacOS; on Ubuntu, it can spit
    # out errors and exceptions about broken pipes, 404 files, etc.

    if sys.platform == "darwin":
      browser_url = http_url
    # Try the default browser.

    if args.browser:
      webbrowser.open(browser_url)
    # Wait for server.
    print()
    print("Serving and waiting ... use Ctrl-C to exit.")
    httpd.serve_forever()


if __name__ == "__main__":
    gen_script = join(dirname(abspath(__file__)), "gen_sphinx.py")
    preview_main(gen_script=gen_script, default_port=8001)
