#!/usr/bin/env python3
"""
Simple HTTP server for Tudat WASM test suite.
Serves files with correct MIME types for WASM.
"""
import http.server
import socketserver
import webbrowser
import os
import sys

PORT = 8832  # TUD on phone keypad (T=8, U=8, D=3) + 2 for "to"
DIRECTORY = os.path.dirname(os.path.abspath(__file__))

class ReusableTCPServer(socketserver.TCPServer):
    """TCP server that allows address reuse for clean restarts."""
    allow_reuse_address = True

class WasmHandler(http.server.SimpleHTTPRequestHandler):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=DIRECTORY, **kwargs)

    def end_headers(self):
        # Required headers for WASM with SharedArrayBuffer support
        # Using 'credentialless' instead of 'require-corp' to allow CDN resources
        self.send_header('Cross-Origin-Opener-Policy', 'same-origin')
        self.send_header('Cross-Origin-Embedder-Policy', 'credentialless')
        super().end_headers()

    def guess_type(self, path):
        if path.endswith('.wasm'):
            return 'application/wasm'
        if path.endswith('.js'):
            return 'application/javascript'
        return super().guess_type(path)

def main():
    os.chdir(DIRECTORY)

    with ReusableTCPServer(("", PORT), WasmHandler) as httpd:
        url = f"http://localhost:{PORT}"
        print(f"\n{'='*60}")
        print(f"  Tudat WASM Embind Test Suite")
        print(f"{'='*60}")
        print(f"\n  Server running at: {url}")
        print(f"  Press Ctrl+C to stop\n")
        print(f"{'='*60}\n")

        # Open browser
        webbrowser.open(url)

        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            print("\n\nServer stopped.")
            sys.exit(0)

if __name__ == "__main__":
    main()
