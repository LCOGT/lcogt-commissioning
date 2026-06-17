"""Simple SciOps text paste processor GUI.

Paste text into the box (Cmd+V). The text is processed immediately,
replaced in the textbox, and copied to the system clipboard.
"""

from __future__ import annotations

import tkinter as tk
from tkinter import ttk

from sciopstools import convert_proposals_to_review_spreadsheet


def process_input_text(input_text: str) -> str:
    """Transform pasted text.

    Edit this function with your actual transformation logic.
    Current behavior normalizes trailing spaces and strips outer blank lines.
    """
    
    return convert_proposals_to_review_spreadsheet.convertproposals(input_text)


class TextPasteProcessorApp:
    def __init__(self, root: tk.Tk) -> None:
        self.root = root
        self.root.title("SciOps Text Paste Processor")
        self.root.geometry("900x500")

        container = ttk.Frame(root, padding=12)
        container.pack(fill=tk.BOTH, expand=True)

        label = ttk.Label(
            container,
            text="Paste text below (Cmd+V). It will be processed and copied automatically.",
        )
        label.pack(anchor=tk.W, pady=(0, 8))

        self.text = tk.Text(container, wrap=tk.WORD, undo=True)
        self.text.pack(fill=tk.BOTH, expand=True)
        self.text.focus_set()

        # Handle keyboard paste and virtual paste event.
        self.text.bind("<Command-v>", self._on_paste_event, add=True)
        self.text.bind("<Control-v>", self._on_paste_event, add=True)
        self.text.bind("<<Paste>>", self._on_paste_event, add=True)

    def _on_paste_event(self, event: tk.Event) -> str | None:
        # Allow Tk to perform the actual paste first, then process.
        self.root.after_idle(self._process_current_text)
        return None

    def _process_current_text(self) -> None:
        raw_text = self.text.get("1.0", tk.END)
        processed_text = process_input_text(raw_text)

        self.text.delete("1.0", tk.END)
        self.text.insert("1.0", processed_text)

        self.root.clipboard_clear()
        self.root.clipboard_append(processed_text)
        self.root.update_idletasks()


def main() -> None:
    root = tk.Tk()
    TextPasteProcessorApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
