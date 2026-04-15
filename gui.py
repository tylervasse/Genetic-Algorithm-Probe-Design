# gui.py — Graphical UI for designing RNA probe structure constraints.
# Supports variable-length alternative sequence blocks and interactive
# dot-bracket structure editing with optional ViennaRNA layout.

import json
import math
import tkinter as tk
from tkinter import font as tkfont
from tkinter import ttk, messagebox, simpledialog, filedialog
import csv
import os
from pathlib import Path
from datetime import datetime


# Optional ViennaRNA for nice coordinates
try:
    import RNA
    HAS_VIENNA = True
except Exception:
    HAS_VIENNA = False

# ---------- Helpers from your pipeline ----------
from dna_utils import rnatoDNA, revcomp

def uneditableProbes(miRNA):
    miDNA = rnatoDNA(miRNA)
    placeholder = revcomp(miDNA) + "TTTT"
    probe = revcomp(placeholder)[:-9]
    probe_end = probe[-7:]
    probe = "AAAA" + revcomp(probe_end) + "TATTAA" + probe
    return probe[17:]

# ---------- Dot-bracket parsing ----------
def dotbracket_pairs(ss: str):
    stack, pairs = [], []
    for i, c in enumerate(ss):
        if c == "(":
            stack.append(i)
        elif c == ")":
            if stack:
                j = stack.pop()
                pairs.append((j, i))
    return pairs

# ---------- Layout ----------
def layout_coords(ss: str, W=900, H=520, margin=50):
    W = max(W, 200)
    H = max(H, 160)
    n = len(ss)
    if n == 0:
        return []
    if HAS_VIENNA:
        try:
            RNA.cvar.rna_plot_type = 1
            coords = RNA.get_xy_coordinates(ss)
            xs = [coords.get(i).X for i in range(n)]
            ys = [coords.get(i).Y for i in range(n)]
            minx, maxx = min(xs), max(xs)
            miny, maxy = min(ys), max(ys)
            spanx = max(1e-6, maxx - minx)
            spany = max(1e-6, maxy - miny)
            sx = max(1e-6, (W - 2*margin) / spanx)
            sy = max(1e-6, (H - 2*margin) / spany)
            s = min(sx, sy)
            X = [margin + (x - minx)*s for x in xs]
            Y = [margin + (y - miny)*s for y in ys]
            return list(zip(X, Y))
        except Exception:
            pass
    # Circular fallback
    r = max(10, min(W, H)/2 - margin)
    cx, cy = W/2, H/2
    return [(cx + r*math.cos(2*math.pi*i/n - math.pi/2),
             cy + r*math.sin(2*math.pi*i/n - math.pi/2)) for i in range(n)]


# ---------- UI ----------
class ConstraintUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.group_font = tkfont.Font(family="Helvetica", size=9)

        self.title("RNA Probe Designer — Structure & Constraint Builder")
        self.geometry("1250x1000")

        

        # Core state
        self.desired_structure = ""
        self.n = 0
        self.pairs = []
        self.allowed = []             # per-base allowed set
        self.weights = []             # per-base mismatch weight 0..9

        # Terminal bracket/export settings
        self.start_penalty = 3   # for "(3)" after the leading [+...]
        self.end_penalty   = 9   # for "(9)" after the trailing [...-]
        self.start_span    = 0   # how many leading digits to include inside [+ ... ]
        self.end_span      = 0   # how many trailing digits to include inside [ ... -]

        # Shift-click group feature kept as before (optional)
        self.groups = []              # {start, end, penalty}
        self.group_start = None

        self._ga_cancel = False 

        self._plus_hit = None
        self._minus_hit = None

        self.pending_span_start = None  # stores -1 (plus), 0..n-1 (base), or n (minus)

        # Entire-section alternative tokens: replaces bases[start..end] with one [ALT1|ALT2] token
        # Each: {"start": int (0..n-1), "end": int (start..n-1), "alts": ["AAA", "AA", ...]}
        self.section_alts = []

        self._build_layout()

        # Seed but DON'T render yet
        example = "....(((((....((((((...........))))))....)))))"
        self.ss_entry.insert(0, example)

        # Ensure the first render uses real canvas size
        self.after_idle(self._first_real_draw)

        #self._render_from_entry()




    # ----- Layout -----
    def _build_layout(self):
        # --- Status bar: always visible at the bottom ---
        self.status = tk.StringVar(value="Ready")
        statusbar = ttk.Frame(self)
        statusbar.pack(side=tk.BOTTOM, fill=tk.X)
        ttk.Label(statusbar, textvariable=self.status, anchor="w").pack(side=tk.LEFT, padx=8, pady=4)

        top = ttk.Frame(self)
        top.pack(side=tk.TOP, fill=tk.X, padx=10, pady=8)

        ttk.Label(top, text="Desired structure (dot-bracket):").pack(side=tk.LEFT)
        self.ss_entry = ttk.Entry(top, width=90)
        self.ss_entry.pack(side=tk.LEFT, padx=8)
        ttk.Button(top, text="Render", command=self._render_from_entry).pack(side=tk.LEFT, padx=4)
        ttk.Button(top, text="Clear", command=self._clear_all).pack(side=tk.LEFT)

        # ===== NEW: use a PanedWindow so the right-hand side is resizable =====
        mid = ttk.PanedWindow(self, orient=tk.HORIZONTAL)
        mid.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=10, pady=8)

        # Left pane: canvas container
        left_frame = ttk.Frame(mid)
        mid.add(left_frame, weight=3)  # weight is for ttk theme support; some themes ignore it

        self.canvas = tk.Canvas(left_frame, bg="white")
        self.canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.canvas.bind("<Button-1>", self._on_canvas_click)

        self.canvas.bind("<Configure>", self._on_canvas_resize)

        # --- Drag/pan state ---
        self._drag_ctx = None       # {"kind": "group"|"section", "idx": int, "items": [canvas_ids], "sx":int, "sy":int}
        self._space_pan = False     # hold Space to pan with left-drag

        # Middle mouse pan
        self.canvas.bind("<ButtonPress-2>", self._scan_mark)
        self.canvas.bind("<B2-Motion>", self._scan_drag)

        # Space-to-pan
        self.bind("<KeyPress-space>", lambda e: setattr(self, "_space_pan", True))
        self.bind("<KeyRelease-space>", lambda e: setattr(self, "_space_pan", False))
        self.canvas.bind("<B1-Motion>", self._maybe_scan_drag)  # only used when space is held

        # Right pane: tools panel (width is initial size; user can drag sash)
        right = ttk.Frame(mid, width=360)
        mid.add(right, weight=1)
        # ===== END NEW LAYOUT BLOCK =====

        self.ims_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(top, text="iMS probe shortcut",
                        variable=self.ims_var,
                        command=self._on_toggle_ims).pack(side=tk.LEFT, padx=12)

        # Per-base tools
        box = ttk.LabelFrame(right, text="Per-base tools")
        box.pack(fill=tk.X, pady=6)
        ttk.Button(box, text="Reset all seq constraints to N", command=self._reset_all_seq).pack(fill=tk.X, pady=2)
        ttk.Button(box, text="Reset all weights to 0", command=self._reset_all_weights).pack(fill=tk.X, pady=2)

        # Section alts (variable-length tokens)
        abox = ttk.LabelFrame(right, text="Section alts")
        abox.pack(fill=tk.X, pady=8)
        ttk.Label(abox, text="Replace bases [start→end] with [ORIGINAL|ALT1|ALT2].").pack(anchor="w")

        ttk.Button(abox, text="Add section alt…", command=self._add_section_alt_dialog).pack(fill=tk.X, pady=3)


        ttk.Label(abox, text="Current section alts: (start→end : [options])").pack(anchor="w", pady=(6, 2))
        self.section_list = tk.Listbox(abox, height=5)
        self.section_list.pack(fill=tk.X, pady=2)
        ttk.Button(abox, text="Clear section alts", command=self._clear_section_alts).pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)
        ttk.Button(abox, text="Delete selected section alt", command=self._delete_selected_section).pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)



        # Gap groups (optional, unchanged)
        gbox = ttk.LabelFrame(right, text="Gap groups (optional)")
        gbox.pack(fill=tk.X, pady=8)
        ttk.Label(gbox, text="Create: Shift+Click start location → Click end location").pack(anchor="w")
        self.group_list = tk.Listbox(gbox, height=5)
        self.group_list.pack(fill=tk.X, pady=2)
        btns = ttk.Frame(gbox)
        btns.pack(fill=tk.X)
        ttk.Button(btns, text="Clear groups", command=self._clear_groups).pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)
        ttk.Button(btns, text="Delete selected", command=self._delete_selected_group).pack(side=tk.LEFT, expand=True, fill=tk.X, padx=2)


        # Export
        ebox = ttk.LabelFrame(right, text="Export")
        ebox.pack(fill=tk.BOTH, expand=True, pady=8)

        self.out_desired = tk.Text(ebox, height=2, wrap="none")
        self.out_seq = tk.Text(ebox, height=3, wrap="none")
        self.out_struct = tk.Text(ebox, height=2, wrap="none")
        self.out_groups = tk.Text(ebox, height=4, wrap="none")

        ttk.Button(ebox, text="Export strings", command=self._export_all).pack(fill=tk.X, pady=4)
        ttk.Label(ebox, text="desired_structure:").pack(anchor="w")
        self.out_desired.pack(fill=tk.X)
        ttk.Label(ebox, text="seq_constraint:").pack(anchor="w")
        self.out_seq.pack(fill=tk.X)
        ttk.Label(ebox, text="struct_constraints_base:").pack(anchor="w")
        self.out_struct.pack(fill=tk.X)
        # --- New dynamic_penalty export field ---
        ttk.Label(ebox, text="dynamic_penalty:").pack(anchor="w")
        self.out_dyn = tk.Text(ebox, height=2, wrap="none")
        self.out_dyn.pack(fill=tk.X)

        self.run_ga_btn = ttk.Button(ebox, text="Run Genetic Algorithm…",
                             command=self._open_ga_dialog, state="disabled")
        self.run_ga_btn.pack(fill=tk.X, pady=(10, 6))

    def _first_real_draw(self):
        # force geometry calculation
        self.update_idletasks()

        # One-shot: re-render again when canvas first reports a real size
        def _once(_e):
            try:
                self.canvas.unbind("<Configure>", once_id)
            except Exception:
                pass
            self._render_from_entry()

        once_id = self.canvas.bind("<Configure>", _once)

        # Safe first draw (in case Configure already fired)
        self._render_from_entry()

    def _on_canvas_resize(self, event):
        """
        Redraw the structure whenever the canvas is resized
        (e.g., when the right-hand pane width changes).
        """
        # If no structure yet, nothing to do
        if not self.desired_structure:
            return

        # Re-draw using the new canvas width/height
        self._draw()


    def _unit_vec(self, x1, y1, x2, y2):
        dx, dy = x2 - x1, y2 - y1
        d = math.hypot(dx, dy) or 1.0
        return dx / d, dy / d

    def _draw_terminal_marker(self, which="+"):
        if self.n == 0: return
        R = 14
        if which == "+":
            x0, y0 = self.coords[0]
            x1, y1 = self.coords[1] if self.n > 1 else (x0 + 1, y0)
            dx, dy = self._unit_vec(x1, y1, x0, y0)
            x, y = x0 + dx * (R + 18), y0 + dy * (R + 18)
            self._draw_ghost_node(x, y, "+")
            # store a small square for hit testing
            self._plus_hit = (x - 12, y - 12, x + 12, y + 12)
        else:
            xn, yn = self.coords[-1]
            xm, ym = self.coords[-2] if self.n > 1 else (xn - 1, yn)
            dx, dy = self._unit_vec(xm, ym, xn, yn)
            x, y = xn + dx * (R + 18), yn + dy * (R + 18)
            self._draw_ghost_node(x, y, "−")
            self._minus_hit = (x - 12, y - 12, x + 12, y + 12)


    def _draw_ghost_node(self, x, y, label):
        # diamond-shaped ghost node (not a base)
        s = 10
        self.canvas.create_polygon(
            x, y - s, x + s, y, x, y + s, x - s, y,
            fill="#e6f0ff", outline="#88aadd", width=1.5
        )
        self.canvas.create_text(x, y, text=label, font=("Helvetica", 10, "bold"), fill="#335")

    def _hit_plus_or_minus(self, event):
        """Return -1 if + was clicked, self.n if - was clicked, else None."""
        if self._plus_hit:
            x1, y1, x2, y2 = self._plus_hit
            if x1 <= event.x <= x2 and y1 <= event.y <= y2:
                return -1
        if self._minus_hit:
            x1, y1, x2, y2 = self._minus_hit
            if x1 <= event.x <= x2 and y1 <= event.y <= y2:
                return self.n
        return None
    
    def _add_section_alt_dialog(self):
        if self.n == 0:
            messagebox.showwarning("No structure", "Render a structure first.")
            return

        s = simpledialog.askinteger("Section start",
                                    f"Start base index (1 - {self.n}):",
                                    parent=self, minvalue=1, maxvalue=self.n)
        if s is None:
            return
        e = simpledialog.askinteger("Section end",
                                    f"End base index (≥ {s}, ≤ {self.n}):",
                                    parent=self, minvalue=s, maxvalue=self.n)
        if e is None:
            return

        opts = simpledialog.askstring("Alternatives",
                                    "Enter options separated by '|', e.g. AAA|AA|TATTAA",
                                    parent=self)
        if not opts:
            return
        alts = [a.strip().upper() for a in opts.split("|") if a.strip()]
        if not alts:
            return

        # Store 0-based inclusive indices
        self.section_alts.append({"start": s-1, "end": e-1, "alts": alts})
        # Keep deterministic order, merge obvious duplicates/normalize
        self.section_alts.sort(key=lambda b: (b["start"], b["end"], "|".join(b["alts"])))
        self._refresh_section_list()
        self._draw()

    def _refresh_section_list(self):
        self.section_list.delete(0, tk.END)
        for b in self.section_alts:
            self.section_list.insert(
                tk.END,
                f"{b['start']+1}→{b['end']+1}  :  [{ '|'.join(b['alts']) }]"
            )

    def _delete_selected_section(self):
        sel = self.section_list.curselection()
        if not sel:
            return
        del self.section_alts[sel[0]]
        self._refresh_section_list()
        self._draw()

    def _clear_section_alts(self):
        if not self.section_alts:
            return
        if messagebox.askyesno("Clear section alts", "Remove all section alts?"):
            self.section_alts = []
            self._refresh_section_list()
            self._draw()

    def _refresh_alt_list(self):
        pass

    def _on_toggle_ims(self):
        # Only act on check; then immediately uncheck so it can be used again.
        if self.ims_var.get():
            self._run_ims_shortcut()
            self.ims_var.set(False)

    def _run_ims_shortcut(self):
        # Ask for miRNA (RNA alphabet)
        miRNA = simpledialog.askstring(
            "iMS probe shortcut",
            "Enter miRNA (RNA, e.g., UAGC...):",
            parent=self
        )
        if not miRNA:
            return

        # remember for GA runs (no prompt in GA dialog)
        self.current_miRNA = miRNA.strip().upper()

        # 1) Desired structure (fixed iMS template, length = 36)
        ss = "....((((((....................))))))"
        self.ss_entry.delete(0, tk.END)
        self.ss_entry.insert(0, ss)
        self._render_from_entry()  # sets self.n, self.pairs, allocates allowed/weights

        # 2) Compute probes
        try:
            u1 = uneditableProbes(miRNA)  # FULL probe (will become ORIGINAL)
            u2 = u1[:-1]                  # shorter alt
        except Exception as e:
            messagebox.showerror("Error", f"Failed to build probes: {e}")
            return

        # 3) Place u1 so that it ends at the last base (no change to structure length)
        L = len(u1)
        if L > self.n:
            messagebox.showerror(
                "Error",
                f"Probe of length {L} is too long for the iMS template of length {self.n}."
            )
            return

        start = self.n - L          # inclusive start index
        end   = self.n - 1          # inclusive end index

        # 4) Write u1 into the per-base constraints
        #    (earlier bases remain 'N' by default)
        for i, ch in enumerate(u1):
            self.allowed[start + i] = {ch}   # lock ORIGINAL to u1

        # 5) Record the section alt with ONLY the shorter probe (u2)
        #    exporter will use ORIGINAL = left-hand bases (u1), so you get [u1|u2]
        #    Remove any overlapping section_alts first
        self.section_alts = [
            s for s in getattr(self, "section_alts", [])
            if not (max(start, s["start"]) <= min(end, s["end"]))
        ]
        self.section_alts.append({"start": start, "end": end, "alts": [u2]})
        self.section_alts.sort(key=lambda b: (b["start"], b["end"], "|".join(b["alts"])))
        self._refresh_section_list()

        # 6) Struct weights to match your reference mid-string
        #    struct_constraints_base target: "33333" + "000000000111111100000000000008" + "9"
        mid = "000000000111111100000000000008"
        full_digits = "3" * 5 + mid + "9"    # length should be 36
        # In case the template is ever changed, trim/pad to n conservatively
        if len(full_digits) != self.n:
            if len(full_digits) > self.n:
                full_digits = full_digits[:self.n]
            else:
                full_digits = full_digits + "0" * (self.n - len(full_digits))
        self.weights = [int(ch) for ch in full_digits]

        # 7) Terminal spans (gap groups) to emulate …
        self.groups = []
        if self.n > 0:
            # + through base 5 (indices 0..4)
            self.groups.append({"start": -1, "end": min(4, self.n - 1), "penalty": 3})
            # last base through −
            self.groups.append({"start": self.n - 1, "end": self.n, "penalty": 9})
        self._refresh_group_list()

        # Redraw & export so the strings show immediately
        self._draw()
        self._export_all()


    def _scan_mark(self, event):
        self.canvas.scan_mark(event.x, event.y)

    def _scan_drag(self, event):
        self.canvas.scan_dragto(event.x, event.y, gain=1)

    def _maybe_scan_drag(self, event):
        # Only pan with left button while space is held; otherwise let normal click/drag handlers do their thing.
        if self._space_pan and self._drag_ctx is None:
            self._scan_drag(event)

    def _drag_start(self, event):
        item = self.canvas.find_withtag("current")
        if not item:
            return
        tags = self.canvas.gettags(item[0])
        kind, idx = None, None
        for t in tags:
            if t.startswith("group_"):
                suf = t[6:]
                if suf.isdigit():
                    kind, idx = "group", int(suf)
                    break
            if t.startswith("section_"):
                suf = t[8:]
                if suf.isdigit():
                    kind, idx = "section", int(suf)
                    break
        if kind is None:
            return
        bundle_tag = f"{kind}_{idx}"
        items = list(self.canvas.find_withtag(bundle_tag))
        self._drag_ctx = {
            "kind": kind, "idx": idx, "items": items,
            "sx": event.x, "sy": event.y, "px": event.x, "py": event.y
        }


    def _drag_move(self, event):
        if not self._drag_ctx: return
        dx = event.x - self._drag_ctx["sx"]
        dy = event.y - self._drag_ctx["sy"]
        for it in self._drag_ctx["items"]:
            self.canvas.move(it, dx, dy)
        self._drag_ctx["sx"], self._drag_ctx["sy"] = event.x, event.y

    def _drag_end(self, event):
        if not self._drag_ctx: return
        total_dx = event.x - self._drag_ctx["px"]
        total_dy = event.y - self._drag_ctx["py"]
        kind = self._drag_ctx["kind"]; idx = self._drag_ctx["idx"]

        if kind == "group":
            g = self.groups[idx]
            vo = g.get("vis_offset", {"dx": 0, "dy": 0})
            vo["dx"] += total_dx; vo["dy"] += total_dy
            g["vis_offset"] = vo
        else:
            b = self.section_alts[idx]
            vo = b.get("vis_offset", {"dx": 0, "dy": 0})
            vo["dx"] += total_dx; vo["dy"] += total_dy
            b["vis_offset"] = vo

        self._drag_ctx = None
        # Redraw once to normalize positions (and clamp if needed)
        self._draw()

    def _open_progress(self, message="Working..."):
        # Reset cancellation flag whenever a new GA run starts
        self._ga_cancel = False

        # Avoid duplicates
        if getattr(self, "_progress_win", None) and self._progress_win.winfo_exists():
            try:
                self._progress_msg.set(message)
            except Exception:
                pass
            return

        # Create, but keep it hidden while we lay it out
        self._progress_win = tk.Toplevel(self)
        self._progress_win.withdraw()  # hide while positioning
        self._progress_win.title("Running Genetic Algorithm")
        self._progress_win.transient(self)
        self._progress_win.grab_set()   # modal
        self._progress_win.resizable(False, False)

        frm = ttk.Frame(self._progress_win, padding=12)
        frm.pack(fill=tk.BOTH, expand=True)

        self._progress_msg = tk.StringVar(value=message)
        ttk.Label(frm, textvariable=self._progress_msg).pack(anchor="w", pady=(0, 8))

        self._progress_bar = ttk.Progressbar(frm, mode="indeterminate", length=320)
        self._progress_bar.pack(fill=tk.X)
        self._progress_bar.start(12)  # smaller = faster animation

        btns = ttk.Frame(frm)
        btns.pack(fill=tk.X, pady=(10, 0))

        # --- Cancel button: request GA cancellation ---
        def on_cancel():
            self._ga_cancel = True
            # Update message so user knows it's doing something
            try:
                self._progress_msg.set("Cancelling… waiting for current GA step to finish.")
            except Exception:
                pass

        ttk.Button(btns, text="Cancel", command=on_cancel).pack(side=tk.RIGHT, padx=6)

        # Optional: a "Hide" button (doesn't cancel GA, just closes dialog)
        ttk.Button(btns, text="Hide", command=self._close_progress).pack(side=tk.RIGHT)

        # --- Center over parent window ---
        self.update_idletasks()
        self._progress_win.update_idletasks()

        parent_x = self.winfo_rootx()
        parent_y = self.winfo_rooty()
        parent_w = self.winfo_width()
        parent_h = self.winfo_height()

        win_w = self._progress_win.winfo_width()
        win_h = self._progress_win.winfo_height()

        if parent_w <= 1 or parent_h <= 1:
            # Fallback: center on screen
            screen_w = self._progress_win.winfo_screenwidth()
            screen_h = self._progress_win.winfo_screenheight()
            x = (screen_w - win_w) // 2
            y = (screen_h - win_h) // 2
        else:
            x = parent_x + (parent_w - win_w) // 2
            y = parent_y + (parent_h - win_h) // 2

        self._progress_win.geometry(f"+{x}+{y}")
        self._progress_win.deiconify()   # show at the centered position
        self._progress_win.lift(self)    # bring above the main window



    def _update_progress(self, message=None):
        if getattr(self, "_progress_win", None) and self._progress_win.winfo_exists():
            if message is not None:
                try:
                    self._progress_msg.set(message)
                except Exception:
                    pass

    def _close_progress(self):
        try:
            if getattr(self, "_progress_bar", None):
                self._progress_bar.stop()
            if getattr(self, "_progress_win", None) and self._progress_win.winfo_exists():
                self._progress_win.grab_release()
                self._progress_win.destroy()
        except Exception:
            pass
        finally:
            self._progress_win = None
            self._progress_bar = None
            self._progress_msg = None


    def _open_ga_dialog(self):
        desired = self.out_desired.get("1.0", "end-1c").strip()
        seq_constraint = self.out_seq.get("1.0", "end-1c").strip()
        dynamic_penalty = self.out_dyn.get("1.0", "end-1c").strip()

        if not desired or not seq_constraint or not dynamic_penalty:
            messagebox.showerror("Missing export",
                                "Please Export strings first so we can pass the current design to the GA.")
            return

        win = tk.Toplevel(self)
        win.title("Run Genetic Algorithm")
        win.transient(self)
        win.grab_set()

        row = 0

        # --- NEW: checkbox to enable/disable Tm constraint ---
        use_tm_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            win,
            text="Use melting temperature constraint",
            variable=use_tm_var
        ).grid(row=row, column=0, columnspan=2, sticky="w", padx=8, pady=(8, 4))
        row += 1

        # Target Tm widgets (default 49°C)
        tmt_var = tk.DoubleVar(value=49.0)
        tm_label = ttk.Label(win, text="Target melting temperature (°C):")
        tm_label.grid(row=row, column=0, sticky="w", padx=8, pady=6)
        tm_spin = ttk.Spinbox(
            win,
            from_=0,
            to=100,
            increment=0.1,
            textvariable=tmt_var,
            width=12
        )
        tm_spin.grid(row=row, column=1, padx=8, pady=6, sticky="w")
        row += 1

        # Enable/disable Tm widgets when box is toggled
        def _toggle_tm_widgets(*_):
            if use_tm_var.get():
                tm_label.state(["!disabled"])
                tm_spin.state(["!disabled"])
            else:
                tm_label.state(["disabled"])
                tm_spin.state(["disabled"])

        use_tm_var.trace_add("write", _toggle_tm_widgets)
        _toggle_tm_widgets()

        # --- NEW: penalize off-target structures (alt structures) toggle ---
        pen_off_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            win,
            text="Penalize off-target structures (extra suboptimal structures)",
            variable=pen_off_var
        ).grid(row=row, column=0, columnspan=2, sticky="w", padx=8, pady=4)
        row += 1

        ttk.Label(win, text="Population size:").grid(row=row, column=0, sticky="w", padx=8, pady=6)
        pop_var = tk.IntVar(value=250)
        ttk.Spinbox(win, from_=10, to=10000, increment=10, textvariable=pop_var, width=12).grid(row=row, column=1, padx=8, pady=6, sticky="w"); row += 1

        ttk.Label(win, text="Generations:").grid(row=row, column=0, sticky="w", padx=8, pady=6)
        gen_var = tk.IntVar(value=100)
        ttk.Spinbox(win, from_=1, to=10000, textvariable=gen_var, width=12).grid(row=row, column=1, padx=8, pady=6, sticky="w"); row += 1

        ttk.Label(win, text="Runs:").grid(row=row, column=0, sticky="w", padx=8, pady=6)
        runs_var = tk.IntVar(value=10)
        ttk.Spinbox(win, from_=1, to=1000, textvariable=runs_var, width=12).grid(row=row, column=1, padx=8, pady=6, sticky="w"); row += 1

        # after the "Runs:" Spinbox
        ttk.Label(win, text="Runs at a time:").grid(row=row, column=0, sticky="w", padx=8, pady=6)
        batch_var = tk.IntVar(value=1)  # default 1
        ttk.Spinbox(win, from_=1, to=32, textvariable=batch_var, width=12).grid(row=row, column=1, padx=8, pady=6, sticky="w"); row += 1


        ttk.Label(win, text="[Na+] (M):").grid(row=row, column=0, sticky="w", padx=8, pady=6)
        na_var = tk.DoubleVar(value=0.30)
        ttk.Spinbox(win, from_=0.0, to=1.0, increment=0.01, textvariable=na_var, width=12).grid(row=row, column=1, padx=8, pady=6, sticky="w"); row += 1

        ttk.Label(win, text="[Mg2+] (M):").grid(row=row, column=0, sticky="w", padx=8, pady=6)
        mg_var = tk.DoubleVar(value=0.0)
        ttk.Spinbox(win, from_=0.0, to=1.0, increment=0.01, textvariable=mg_var, width=12).grid(row=row, column=1, padx=8, pady=6, sticky="w"); row += 1

        # Placeholder options (all gated by generate-placeholder checkbox)
        ph_frame = ttk.LabelFrame(win, text="Placeholder strand options")
        ph_frame.grid(row=row, column=0, columnspan=2, sticky="ew", padx=8, pady=(8, 6)); row += 1

        gen_ph_var = tk.BooleanVar(value=True)
        auto_ph_var = tk.BooleanVar(value=True)
        ph_tm_var = tk.DoubleVar(value=55.0)

        gen_chk  = ttk.Checkbutton(ph_frame, text="Generate placeholder strand", variable=gen_ph_var)
        auto_chk = ttk.Checkbutton(ph_frame, text="Auto-calc placeholder Tm", variable=auto_ph_var)
        ph_lbl   = ttk.Label(ph_frame, text="Desired placeholder Tm (°C):")
        ph_spin  = ttk.Spinbox(ph_frame, from_=0, to=100, increment=0.1, textvariable=ph_tm_var, width=12)

        gen_chk.grid(row=0, column=0, sticky="w", padx=8, pady=4)
        auto_chk.grid(row=1, column=0, sticky="w", padx=8, pady=4)
        ph_lbl.grid(row=2, column=0, sticky="w", padx=8, pady=4)
        ph_spin.grid(row=2, column=1, sticky="w", padx=8, pady=4)

        def _apply_ph_enable_states(*_):
            # If generate-placeholder is OFF: disable auto checkbox and the spinbox
            if not gen_ph_var.get():
                auto_chk.state(["disabled"])
                ph_spin.state(["disabled"])
            else:
                # Generate-placeholder ON: enable auto checkbox; spinbox depends on auto flag
                auto_chk.state(["!disabled"])
                if auto_ph_var.get():
                    ph_spin.state(["disabled"])
                else:
                    ph_spin.state(["!disabled"])

        gen_ph_var.trace_add("write", _apply_ph_enable_states)
        auto_ph_var.trace_add("write", _apply_ph_enable_states)
        _apply_ph_enable_states()

        # ===== NEW: CSV output options =====
        csv_frame = ttk.LabelFrame(win, text="CSV output")
        csv_frame.grid(row=row, column=0, columnspan=2, sticky="ew", padx=8, pady=(8, 6))
        row += 1

        # Default folder: Downloads
        default_dir = str(self._get_downloads_dir())
        csv_dir_var = tk.StringVar(value=default_dir)
        csv_name_var = tk.StringVar(value="ga_results.csv")

        ttk.Label(csv_frame, text="File name:").grid(row=0, column=0, sticky="w", padx=8, pady=4)
        ttk.Entry(csv_frame, textvariable=csv_name_var, width=28).grid(row=0, column=1, sticky="w", padx=8, pady=4)

        ttk.Label(csv_frame, text="Save folder:").grid(row=1, column=0, sticky="w", padx=8, pady=4)
        ttk.Entry(csv_frame, textvariable=csv_dir_var, width=28).grid(row=1, column=1, sticky="w", padx=8, pady=4)

        def browse_folder():
            folder = filedialog.askdirectory(parent=win, initialdir=csv_dir_var.get() or default_dir)
            if folder:
                csv_dir_var.set(folder)

        ttk.Button(csv_frame, text="Browse…", command=browse_folder).grid(row=1, column=2, sticky="w", padx=6, pady=4)


        # Run & Cancel
        btns = ttk.Frame(win); btns.grid(row=row, column=0, columnspan=2, sticky="e", padx=8, pady=8)

        def do_run():
            # ----- NEW: validate CSV name + folder -----
            name = (csv_name_var.get() or "").strip()
            folder = (csv_dir_var.get() or "").strip()

            if not name:
                messagebox.showerror("Missing CSV name", "Please enter a CSV file name (e.g., results.csv).")
                return
            if not folder:
                messagebox.showerror("Missing folder", "Please choose where to save the CSV.")
                return

            # Ensure .csv extension
            if not name.lower().endswith(".csv"):
                name += ".csv"

            out_csv_path = Path(folder) / name


            if gen_ph_var.get() and not getattr(self, "current_miRNA", ""):
                messagebox.showerror(
                    "Missing miRNA",
                    "To generate a placeholder strand you must set the miRNA first via the iMS probe shortcut."
                )
                return

            # Show progress dialog
            self._open_progress("Running genetic algorithm…")

            import threading
            threading.Thread(
                target=self._run_ga_in_thread,
                args=(
                    bool(use_tm_var.get()),
                    bool(pen_off_var.get()),
                    float(tmt_var.get()),
                    desired,
                    seq_constraint,
                    dynamic_penalty,
                    int(pop_var.get()),
                    int(gen_var.get()),
                    int(runs_var.get()),
                    float(na_var.get()),
                    float(mg_var.get()),
                    bool(gen_ph_var.get()),
                    bool(auto_ph_var.get()),
                    float(ph_tm_var.get()),
                    int(batch_var.get()),
                    out_csv_path,  # NEW (last)
                ),
                daemon=True
            ).start()


            win.destroy()


        ttk.Button(btns, text="Run", command=do_run).pack(side=tk.RIGHT, padx=6)
        ttk.Button(btns, text="Cancel", command=win.destroy).pack(side=tk.RIGHT)

        # >>> ADD THIS BLOCK <<<
        self.update_idletasks()
        win.update_idletasks()

        parent_x = self.winfo_rootx()
        parent_y = self.winfo_rooty()
        parent_w = self.winfo_width()
        parent_h = self.winfo_height()
        win_w = win.winfo_width()
        win_h = win.winfo_height()

        x = parent_x + (parent_w - win_w) // 2
        y = parent_y + (parent_h - win_h) // 2
        win.geometry(f"+{x}+{y}")


    def _run_ga_in_thread(self,
                          use_tm_constraint,
                          penalize_off_targets,
                          target_tm,
                          desired,
                          seq_constraint,
                          dynamic_penalty,
                          pop,
                          gens,
                          runs,
                          na,
                          mg,
                          gen_placeholder,
                          auto_ph,
                          ph_tm,
                          batch_size,
                          out_csv_path):   # NEW



        try:
            import genetic_algorithm as GAmod
            mirna = getattr(self, "current_miRNA", "")

            kwargs = dict(
                miRNA=mirna,
                target_melting_temp=target_tm,
                desired_structure=desired,
                seq_constraint=seq_constraint,
                struct_constraints=dynamic_penalty,
                population_size=pop,
                generations=gens,
                runs=runs,
                salt_conc=na,
                mag_conc=mg,
                automatic_ph_mt=auto_ph,
                desired_ph_mt=ph_tm,
                generate_placeholder=gen_placeholder,
                batch_size=batch_size,
                use_tm_constraint=use_tm_constraint,      # existing
                penalize_off_targets=penalize_off_targets # NEW
            )



            # This closure reads the flag that the Cancel button sets
            def cancel_cb():
                return getattr(self, "_ga_cancel", False)

            # Now GA can cooperate with cancellation
            reports = GAmod.GA(cancel_cb=cancel_cb, **kwargs)

            # ----- NEW: write GA results to the user-chosen CSV path -----
            try:
                ga_params = dict(kwargs)  # copy for logging
                ga_params["penalize_off_targets"] = penalize_off_targets
                ga_params["use_tm_constraint"] = use_tm_constraint
                self._write_ga_results_csv(Path(out_csv_path), reports, ga_params)
            except Exception as e:
                err = str(e)
                def show_err(err=err):
                    self._close_progress()
                    messagebox.showerror("GA error", err)
                self.after(0, show_err)

            def open_results():
                self._close_progress()
                if getattr(self, "_ga_cancel", False):
                    # User clicked Cancel: don't open results window
                    messagebox.showinfo("GA cancelled", "Genetic algorithm run was cancelled.")
                    return

                if not reports:
                    messagebox.showinfo("GA finished", "GA completed, but no candidate results were returned.")
                    return
                try:
                    self.status.set(f"GA finished. Results saved to: {out_csv_path}")
                except Exception:
                    pass

                self._open_ga_results_window(reports)

            self.after(0, open_results)

        except Exception as e:
            err = str(e)
            def show_err(err=err):
                self._close_progress()
                messagebox.showerror("GA error", err)
            self.after(0, show_err)

    def _open_ga_results_window(self, reports):
        """Create a paginated window to show one candidate per page with ViennaRNA rendering."""
        self.ga_reports = reports
        self.ga_report_idx = 0

        win = tk.Toplevel(self)
        win.title("Genetic Algorithm Results")
        win.geometry("1250x1000")
        win.transient(self)

        # --- Top: navigation ---
        nav = ttk.Frame(win); nav.pack(side=tk.TOP, fill=tk.X, padx=8, pady=6)
        self._ga_idx_var = tk.IntVar(value=1)
        self._ga_total   = len(reports)

        ttk.Button(nav, text="◀ Prev", command=lambda: self._ga_nav(-1, win)).pack(side=tk.LEFT, padx=4)
        ttk.Label(nav, textvariable=self._ga_idx_var).pack(side=tk.LEFT)
        ttk.Label(nav, text=f"/ {self._ga_total}").pack(side=tk.LEFT)
        ttk.Button(nav, text="Next ▶", command=lambda: self._ga_nav(+1, win)).pack(side=tk.LEFT, padx=4)

        # --- Middle: canvas for structure ---
        self.ga_canvas = tk.Canvas(win, bg="white", height=450)
        self.ga_canvas.pack(side=tk.TOP, fill=tk.X, padx=8, pady=6)

        # Ensure first render uses real widget size
        win.update_idletasks()

        # >>> ADD THIS BLOCK <<<
        self.update_idletasks()
        parent_x = self.winfo_rootx()
        parent_y = self.winfo_rooty()
        parent_w = self.winfo_width()
        parent_h = self.winfo_height()
        win_w = win.winfo_width()
        win_h = win.winfo_height()

        x = parent_x + (parent_w - win_w) // 2
        y = parent_y + (parent_h - win_h) // 2
        win.geometry(f"+{x}+{y}")
        # <<< END BLOCK >>>

        self.ga_canvas.after_idle(lambda: self._render_ga_result_page(win))

        def _ga_once(_e):
            try:
                self.ga_canvas.unbind("<Configure>", once_id)
            except Exception:
                pass
            self._render_ga_result_page(win)

        once_id = self.ga_canvas.bind("<Configure>", _ga_once)

        # --- Bottom: details ---
        info = ttk.Frame(win); info.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=8, pady=6)

        self._ga_seq_txt  = tk.Text(info, height=4, wrap="word")
        self._ga_seq_txt.pack(fill=tk.X, pady=(0,6))

        grid = ttk.Frame(info); grid.pack(fill=tk.X)
        self._ga_tm_var   = tk.StringVar()
        self._ga_fit_var  = tk.StringVar()
        self._ga_len_var  = tk.StringVar()
        self._ga_alt_var  = tk.StringVar()
        self._ga_ph_var   = tk.StringVar()
        self._ga_phtm_var = tk.StringVar()
        self._ga_mtm_var  = tk.StringVar()

        def row(label, var, r):
            ttk.Label(grid, text=label).grid(row=r, column=0, sticky="w", padx=4, pady=3)
            ttk.Label(grid, textvariable=var).grid(row=r, column=1, sticky="w", padx=4, pady=3)

        # Always-visible rows
        row("Melting Temp (°C):",      self._ga_tm_var,  0)
        row("Fitness:",                self._ga_fit_var, 1)
        row("Length:",                 self._ga_len_var, 2)
        row("Alt structures (count):", self._ga_alt_var, 3)

        # Placeholder-related rows (created as separate label+value widgets so we can hide/show them)
        self._ga_ph_label    = ttk.Label(grid, text="Placeholder:")
        self._ga_ph_value    = ttk.Label(grid, textvariable=self._ga_ph_var)
        self._ga_phtm_label  = ttk.Label(grid, text="Placeholder Tm (°C):")
        self._ga_phtm_value  = ttk.Label(grid, textvariable=self._ga_phtm_var)
        self._ga_mtm_label   = ttk.Label(grid, text="miRNA Tm (°C):")
        self._ga_mtm_value   = ttk.Label(grid, textvariable=self._ga_mtm_var)

        # Grid them in fixed rows; your renderer will grid_remove() if not applicable
        self._ga_ph_label.grid(row=4, column=0, sticky="w", padx=4, pady=3)
        self._ga_ph_value.grid(row=4, column=1, sticky="w", padx=4, pady=3)
        self._ga_phtm_label.grid(row=5, column=0, sticky="w", padx=4, pady=3)
        self._ga_phtm_value.grid(row=5, column=1, sticky="w", padx=4, pady=3)
        self._ga_mtm_label.grid(row=6, column=0, sticky="w", padx=4, pady=3)
        self._ga_mtm_value.grid(row=6, column=1, sticky="w", padx=4, pady=3)

        # Initial render (the after_idle/Configure path above will also trigger, but this is harmless)
        self._render_ga_result_page(win)


    def _ga_nav(self, delta, win):
        new_idx = (self.ga_report_idx + delta) % len(self.ga_reports)
        self.ga_report_idx = new_idx
        self._ga_idx_var.set(new_idx + 1)
        self._render_ga_result_page(win)

    def _render_ga_result_page(self, win):
        rep = self.ga_reports[self.ga_report_idx]

        # --- Melting temp display (handle None) ---
        tm = rep.get("melting_temp")
        tm_str = "N/A" if tm is None else f"{tm:.2f}"
        self._ga_tm_var.set(tm_str)

        seq = rep["sequence"]
        ss  = rep["structure"]

        # Draw structure
        W = self.ga_canvas.winfo_width() or 940
        H = self.ga_canvas.winfo_height() or 450
        self.ga_canvas.delete("all")

        coords = layout_coords(ss, W=W-40, H=H-30, margin=40)
        pairs = dotbracket_pairs(ss)

        # base pairs (blue)
        for i, j in pairs:
            x1, y1 = coords[i]; x2, y2 = coords[j]
            self.ga_canvas.create_line(x1, y1, x2, y2, fill="#69b", width=2)
        # backbone (light gray)
        for i in range(len(seq)-1):
            x1, y1 = coords[i]; x2, y2 = coords[i+1]
            self.ga_canvas.create_line(x1, y1, x2, y2, fill="#aaa", width=1)
        # nodes
        R = 10
        for i, (x, y) in enumerate(coords):
            self.ga_canvas.create_oval(x-R, y-R, x+R, y+R, fill="#eef", outline="#557")
            base = seq[i] if i < len(seq) else "N"
            self.ga_canvas.create_text(x, y-0.5, text=base, font=("Helvetica", 8, "bold"), fill="#223")

        # Fill text fields
        self._ga_seq_txt.delete("1.0", "end")
        self._ga_seq_txt.insert("1.0", f"Sequence (rank {rep['rank']}):\n{seq}")

        self._ga_fit_var.set(f"{rep['fitness']:.2f}")
        self._ga_len_var.set(str(rep['length']))
        self._ga_alt_var.set(str(max(0, rep.get('alt_structures', 0) - 1)))

        # Optional placeholder fields
        has_placeholder = (
            "placeholder" in rep
            and rep.get("placeholder")
            and getattr(self, "current_miRNA", "")
        )

        for widget in (self._ga_ph_label, self._ga_phtm_label, self._ga_mtm_label,
                       self._ga_ph_value, self._ga_phtm_value, self._ga_mtm_value):
            if has_placeholder:
                widget.grid()
            else:
                widget.grid_remove()

        if has_placeholder:
            self._ga_ph_var.set(rep["placeholder"])
            phtm = rep.get("placeholder_tm")
            mtm  = rep.get("mirna_tm")
            self._ga_phtm_var.set("N/A" if phtm is None else f"{phtm:.2f}")
            self._ga_mtm_var.set("N/A" if mtm is None else f"{mtm:.2f}")
        else:
            self._ga_ph_var.set("")
            self._ga_phtm_var.set("")
            self._ga_mtm_var.set("")




    # ----- Render & validation -----
    def _is_valid_dotbracket(self, ss: str) -> bool:
        if any(c not in ".()" for c in ss):
            return False
        bal = 0
        for c in ss:
            if c == "(":
                bal += 1
            elif c == ")":
                bal -= 1
            if bal < 0:
                return False
        return bal == 0

    def _render_from_entry(self):
        ss = self.ss_entry.get().strip()
        if not self._is_valid_dotbracket(ss):
            messagebox.showerror("Invalid dot-bracket", "Unbalanced or invalid characters in structure.")
            return
        self.desired_structure = ss
        self.n = len(ss)
        self.pairs = dotbracket_pairs(ss)

        # Resize state
        if len(self.allowed) != self.n:
            self.allowed = [set("ACGT") for _ in range(self.n)]
        if len(self.weights) != self.n:
            self.weights = [0 for _ in range(self.n)]

        # Clip alt blocks to new range 0..n
        self._draw()

    def _draw(self):
        self.canvas.delete("all")
        W = self.canvas.winfo_width() or 900
        H = self.canvas.winfo_height() or 520

        coords = layout_coords(self.desired_structure, W=W, H=H, margin=60)
        self.coords = coords

        # base pairs
        for i, j in self.pairs:
            x1, y1 = coords[i]; x2, y2 = coords[j]
            self.canvas.create_line(x1, y1, x2, y2, fill="#69b", width=2)

        # backbone
        for i in range(self.n - 1):
            x1, y1 = coords[i]; x2, y2 = coords[i+1]
            self.canvas.create_line(x1, y1, x2, y2, fill="#aaa", width=1)

        # nodes
        R = 14
        for i, (x, y) in enumerate(coords):
            circ = self.canvas.create_oval(x - R, y - R, x + R, y + R, fill="#eef", outline="#557")
            label = self._label_for_allowed(self.allowed[i])
            self.canvas.create_text(x, y-0.5, text=label, font=("Helvetica", 10, "bold"))
            self.canvas.create_text(x, y - (R + 10), text=str(i+1), font=("Helvetica", 8), fill="#666")
            self.canvas.create_text(x, y + (R + 10), text=str(int(self.weights[i])), font=("Helvetica", 9), fill="#a55")


        # groups underlay (optional)
        self._draw_groups_underlay()

        # NEW: draw alt-block badges near insertion sites
        self._draw_alt_badges()

        self.status.set(f"Rendered {self.n} bases. Click a base to edit. Section alts shown as orange pills.")

        self._draw_terminal_marker("+")
        self._draw_terminal_marker("-")


    def _draw_groups_underlay(self):
        if not getattr(self, "coords", None) or not self.groups:
            return

        R = 14            # node radius
        h = 10            # half-height of the band
        pad = 8           # horizontal padding inside the band
        ch = self.canvas.winfo_height() or 520
        cw = self.canvas.winfo_width() or 900

        def label_pos(p):
            if p == -1: return "+"
            if p == self.n: return "−"
            return str(p + 1)

        for i, g in enumerate(self.groups):
            s, e, pen = g["start"], g["end"], int(g["penalty"])

            # Map terminals to geometry indices for coordinates
            gi = 0 if s == -1 else s
            gj = (self.n - 1) if e == self.n else e
            if gi > gj:
                gi, gj = gj, gi

            x1, _ = self.coords[gi]
            x2, _ = self.coords[gj]
            if x2 < x1:
                x1, x2 = x2, x1

            # Compute a preferred vertical placement just outside the nodes (no overlap)
            ys = [self.coords[k][1] for k in range(gi, gj + 1)]
            midx = (x1 + x2) / 2
            above_y = min(ys) - (R + 6)   # a bit higher than before to avoid covering caps
            below_y = max(ys) + (R + 6)
            band_y = above_y if above_y > 10 else below_y
            band_y = max(10, min(ch - 10, band_y))

            # Label & width calc
            label = f"{label_pos(s)}→{label_pos(e)} (gap penalty={pen})"
            text_w = self.group_font.measure(label)

            # Band width must fit the geometric span AND the label (+ padding)
            span_w = abs(x2 - x1)
            need_w = max(span_w, text_w + 2 * pad)

            # Center the band on the span mid; clamp to canvas, while still covering the span
            bx1 = midx - need_w / 2
            bx2 = midx + need_w / 2
            if bx1 < 6:
                shift = 6 - bx1
                bx1 += shift; bx2 += shift
            if bx2 > cw - 6:
                shift = bx2 - (cw - 6)
                bx1 -= shift; bx2 -= shift

            # Ensure the band still covers the geometric x1..x2 of the span
            bx1 = min(bx1, x1)
            bx2 = max(bx2, x2)

            # Colors (optional: different hue for terminal spans)
            is_terminal = (s == -1 or e == self.n)
            fill = "#d7f5b5" if is_terminal else "#f5d7b5"
            outline = "#a6e48a" if is_terminal else "#eeaaaa"

            # --- Apply any saved visual offset (for dragging) ---
            vo = g.get("vis_offset", {"dx": 0, "dy": 0})
            dxo = vo.get("dx", 0); dyo = vo.get("dy", 0)
            bx1o, bx2o, byo = bx1 + dxo, bx2 + dxo, band_y + dyo
            midxo = midx + dxo

            # Draw band + label
            rid = self.canvas.create_rectangle(bx1o, byo - h, bx2o, byo + h, fill=fill, outline=outline)
            tid = self.canvas.create_text(midxo, byo, text=label, font=self.group_font, fill="#333")

            # Tag & bind for drag
            tag = f"group_{i}"
            self.canvas.addtag_withtag("group_tag", rid)
            self.canvas.addtag_withtag("group_tag", tid)
            self.canvas.addtag_withtag(tag, rid)
            self.canvas.addtag_withtag(tag, tid)
            # (Assumes you implemented _drag_start/_drag_move/_drag_end earlier)
            self.canvas.tag_bind(tag, "<ButtonPress-1>", self._drag_start)
            self.canvas.tag_bind(tag, "<B1-Motion>", self._drag_move)
            self.canvas.tag_bind(tag, "<ButtonRelease-1>", self._drag_end)


    def _draw_alt_badges(self):
        if self.n == 0:
            return

        R = 14
        h = 8
        for i, b in enumerate(self.section_alts):
            s, e = max(0, b["start"]), min(self.n - 1, b["end"])
            if s > e:
                continue

            x1, _ = self.coords[s]
            x2, _ = self.coords[e]
            if x2 < x1:
                x1, x2 = x2, x1

            ys = [self.coords[k][1] for k in range(s, e + 1)]
            band_y = min(ys) - (R + 6)

            # Apply any saved visual offset from dragging
            vo = b.get("vis_offset", {"dx": 0, "dy": 0})
            dxo = vo.get("dx", 0); dyo = vo.get("dy", 0)
            x1o, x2o, yo = x1 + dxo, x2 + dxo, band_y + dyo
            midx = (x1o + x2o) / 2

            # Draw band
            rid = self.canvas.create_rectangle(
                x1o, yo - h, x2o, yo + h,
                fill="#ffe7c9", outline="#f3b878"
            )

            # Label like "20->35 (section alt 1)"
            label = f"{s+1}→{e+1} (section alt {i+1})"
            tid = self.canvas.create_text(
                midx, yo,
                text=label, font=("Helvetica", 9, "bold"), fill="#633"
            )

            # Tag & bind for dragging as a bundle
            tag = f"section_{i}"
            self.canvas.addtag_withtag(tag, rid)
            self.canvas.addtag_withtag(tag, tid)
            self.canvas.tag_bind(tag, "<ButtonPress-1>", self._drag_start)
            self.canvas.tag_bind(tag, "<B1-Motion>", self._drag_move)
            self.canvas.tag_bind(tag, "<ButtonRelease-1>", self._drag_end)

    def _pill(self, x, y, text="", fill="#ffd28a", outline="#c98b3a"):
        w, h, r = 26, 14, 7
        o1 = self.canvas.create_oval(x-w/2, y-h/2, x-w/2+2*r, y+h/2, fill=fill, outline=outline)
        o2 = self.canvas.create_oval(x+w/2-2*r, y-h/2, x+w/2, y+h/2, fill=fill, outline=outline)
        r1 = self.canvas.create_rectangle(x-w/2+r, y-h/2, x+w/2-r, y+h/2, fill=fill, outline=outline)
        t1 = self.canvas.create_text(x, y, text=text, font=("Helvetica", 9, "bold"), fill="#633")
        # group the 4 shape ids under a temporary tag returned by caller? Easiest:
        return (o1, o2, r1, t1)


    def _ask_span_penalty(self, which, initial):
        """
        Ask for a penalty (0..9) for the terminal span just set.
        which: "+" or "-"  (only used for dialog label)
        """
        pen = simpledialog.askinteger(
            f"{which} span penalty",
            f"Penalty for {which} span (0–9):",
            parent=self, minvalue=0, maxvalue=9, initialvalue=int(initial)
        )
        return pen if pen is not None else initial


    # ----- Interaction -----
    def _on_canvas_click(self, event):
        # Space+left drag -> panning
        if self._space_pan:
            self._scan_mark(event)
            return
        # Helper: map click to sentinel index: -1 (plus), 0..n-1 (base), n (minus), or None
        def pick_target(event):
            hit = self._hit_plus_or_minus(event)
            if hit is not None:
                return hit
            return self._nearest_node_index(event)

        # Shift-click => start span (can start at +, base, or -)
        if event.state & 0x0001:
            start = pick_target(event)
            if start is None:
                return
            self.pending_span_start = start
            label = "+" if start == -1 else ("−" if start == self.n else f"base {start+1}")
            self.status.set(f"Span start set at {label}. Click the end (base or terminal).")
            return

        # If a span is in progress, finish it with a normal click
        if self.pending_span_start is not None:
            end = pick_target(event)
            if end is None:
                return
            start = self.pending_span_start
            self.pending_span_start = None

            # Normalize (start, end) to ordered pair in terms of *base indices*,
            # while preserving terminal sentinels.
            s, e = start, end
            # If both are bases, order them
            if isinstance(s, int) and isinstance(e, int) and 0 <= s < self.n and 0 <= e < self.n and s > e:
                s, e = e, s

            # Ask once for penalty
            pen = simpledialog.askinteger(
                "Span penalty",
                "Penalty for this span (0–9):",
                parent=self, minvalue=0, maxvalue=9, initialvalue=5
            )
            if pen is None:
                self._draw()
                return

            # Append group
            # Allowed values: s in {-1, 0..n-1}, e in {0..n-1, n}
            # A terminal span is encoded by s=-1 (start terminal) or e=n (end terminal).
            # A mid-span has 0 <= s <= e <= n-1.
            self.groups.append({"start": s, "end": e, "penalty": int(pen)})
            self._refresh_group_list()
            self._draw()
            return

        # Otherwise, normal click edits a base (if any)
        idx = self._nearest_node_index(event)
        if idx is None:
            return
        self._open_base_editor(idx)



    def _nearest_node_index(self, event):
        x, y = event.x, event.y
        best_i, best_d2 = None, 1e9
        for i, (xi, yi) in enumerate(getattr(self, "coords", [])):
            d2 = (xi-x)**2 + (yi-y)**2
            if d2 < best_d2:
                best_d2, best_i = d2, i
        return best_i if best_d2 <= (18**2) else None

    def _open_base_editor(self, idx):
        win = tk.Toplevel(self)
        win.title(f"Base {idx+1} constraints")
        win.transient(self)
        win.grab_set()

        ttk.Label(win, text=f"Base index: {idx+1}").pack(anchor="w", padx=10, pady=(10, 0))

        fr = ttk.LabelFrame(win, text="Allowed bases")
        fr.pack(fill=tk.X, padx=10, pady=6)
        var = {b: tk.BooleanVar(value=(b in self.allowed[idx])) for b in "ACGT"}
        for col, b in enumerate("ACGT"):
            ttk.Checkbutton(fr, text=b, variable=var[b]).grid(row=0, column=col, padx=6, pady=4)

        wfr = ttk.LabelFrame(win, text="Positional mismatch weight (0-9)")
        wfr.pack(fill=tk.X, padx=10, pady=6)
        wvar = tk.IntVar(value=int(self.weights[idx]))
        ttk.Spinbox(wfr, from_=0, to=9, textvariable=wvar, width=4).pack(padx=8, pady=6)

        def apply_and_close():
            new_allowed = {b for b in "ACGT" if var[b].get()}
            if not new_allowed:
                messagebox.showerror("Invalid", "At least one base must be allowed.")
                return
            self.allowed[idx] = new_allowed
            self.weights[idx] = int(wvar.get())
            win.destroy()
            self._draw()

        btns = ttk.Frame(win)
        btns.pack(fill=tk.X, padx=10, pady=8)
        ttk.Button(btns, text="Apply", command=apply_and_close).pack(side=tk.RIGHT, padx=4)
        ttk.Button(btns, text="Cancel", command=win.destroy).pack(side=tk.RIGHT)

        # ----- NEW: center the dialog over the main window -----
        self.update_idletasks()
        win.update_idletasks()
        parent_x = self.winfo_rootx()
        parent_y = self.winfo_rooty()
        parent_w = self.winfo_width()
        parent_h = self.winfo_height()
        win_w = win.winfo_width()
        win_h = win.winfo_height()

        x = parent_x + (parent_w - win_w) // 2
        y = parent_y + (parent_h - win_h) // 2
        win.geometry(f"+{x}+{y}")


    # ----- Utilities -----
    def _label_for_allowed(self, s):
        if s == set("ACGT"): return "N"
        if len(s) == 1: return next(iter(s))
        lab = "".join(sorted(s))
        return lab if len(lab) <= 2 else f"{len(s)}"

    def _reset_all_seq(self):
        if self.n == 0: return
        self.allowed = [set("ACGT") for _ in range(self.n)]
        self._draw()

    def _reset_all_weights(self):
        if self.n == 0: return
        self.weights = [0 for _ in range(self.n)]
        self._draw()

    def _clear_groups(self):
        self.groups = []
        self._refresh_group_list()
        self._draw()

    def _delete_selected_group(self):
        sel = self.group_list.curselection()
        if not sel: return
        del self.groups[sel[0]]
        self._refresh_group_list()
        self._draw()

    def _refresh_group_list(self):
        self.group_list.delete(0, tk.END)
        for g in self.groups:
            s = g["start"]; e = g["end"]; pen = g["penalty"]
            def label_pos(p):
                if p == -1: return "+"
                if p == self.n: return "−"
                return str(p + 1)
            self.group_list.insert(tk.END, f"{label_pos(s)}→{label_pos(e)} : gap penalty={pen}")


    def _clear_all(self):
        self.ss_entry.delete(0, tk.END)
        self.desired_structure = ""
        self.n = 0
        self.pairs = []
        self.allowed = []
        self.weights = []
        self.groups = []
        self.section_alts = []
        self.canvas.delete("all")
        self._refresh_group_list()
        self._refresh_section_list()
        self.status.set("Cleared.")


    def _build_dynamic_penalty(self, struct_digits):
        """
        Build a dynamic_penalty string from self.groups.
        Rules:
        • Each group produces an interval over base indices [s..e], where
            s = max(0, start), e = min(n-1, end), mapping terminals as:
            start == -1 -> s = 0, with prefix_plus=True
            end   ==  n -> e = n-1, with suffix_minus=True
        • Overlaps are merged; if penalties differ, keep the MAX penalty and
            carry forward terminal flags if any interval in the merge had them.
        • Output is a concatenation of plain digits (outside) and bracketed
            chunks (inside), with (+) and (−) markers inside brackets for terminal spans.
        """
        n = self.n
        if n == 0:
            return ""

        # If no groups, return plain digits
        if not self.groups:
            return struct_digits

        # 1) Normalize groups into intervals
        ivals = []
        for g in self.groups:
            gs = g["start"]; ge = g["end"]; pen = int(g["penalty"])
            if gs == -1:
                s = 0
                prefix_plus = True
            else:
                s = int(gs)
                prefix_plus = False
            if ge == n:
                e = n - 1
                suffix_minus = True
            else:
                e = int(ge)
                suffix_minus = False
            if s > e:
                s, e = e, s
            ivals.append([s, e, pen, prefix_plus, suffix_minus])

        # 2) Merge overlaps by start/end, maximizing penalty and OR-ing terminal flags
        ivals.sort(key=lambda x: (x[0], x[1]))
        merged = []
        for s, e, pen, pplus, sminus in ivals:
            if not merged or s > merged[-1][1] + 1:
                merged.append([s, e, pen, pplus, sminus])
            else:
                # overlap/adjacent: merge into last
                ms, me, mpen, mpplus, msminus = merged[-1]
                merged[-1][1] = max(me, e)
                merged[-1][2] = max(mpen, pen)
                merged[-1][3] = (mpplus or pplus)
                merged[-1][4] = (msminus or sminus)

        # 3) Emit string
        out = []
        cursor = 0
        for s, e, pen, pplus, sminus in merged:
            # plain digits before interval
            if cursor < s:
                out.append(struct_digits[cursor:s])

            # bracketed interval
            chunk = struct_digits[s:e+1]
            left = "[+" if pplus else "["
            right_mark = "-" if sminus else ""
            out.append(f"{left}{chunk}{right_mark}]({pen})")

            cursor = e + 1

        # trailing plain digits
        if cursor < n:
            out.append(struct_digits[cursor:])

        return "".join(out)
    
    def _get_downloads_dir(self) -> Path:
        """
        Return the user's Downloads directory in a cross-platform way.
        Works for Windows/macOS/Linux in typical setups.
        """
        home = Path.home()
        downloads = home / "Downloads"
        # Fallback if Downloads doesn't exist for some reason
        return downloads if downloads.exists() else home

    def _write_export_csv(self, desired, seq_str, struct_digits, dynamic_penalty):
        """
        Write a single-row CSV of the final exported results to Downloads.
        Creates a timestamped filename to avoid overwriting.
        """
        downloads_dir = self._get_downloads_dir()
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_path = downloads_dir / f"rna_constraints_export_{ts}.csv"

        # Optional: include extra context that can be useful later
        section_alts_json = json.dumps(self.section_alts, ensure_ascii=False)
        groups_json = json.dumps(self.groups, ensure_ascii=False)

        fieldnames = [
            "timestamp",
            "desired_structure",
            "seq_constraint",
            "struct_constraints_base",
            "dynamic_penalty",
            "n_bases",
            "num_pairs",
            "section_alts_json",
            "groups_json",
        ]

        row = {
            "timestamp": datetime.now().isoformat(timespec="seconds"),
            "desired_structure": desired,
            "seq_constraint": seq_str,
            "struct_constraints_base": struct_digits,
            "dynamic_penalty": dynamic_penalty,
            "n_bases": self.n,
            "num_pairs": len(self.pairs),
            "section_alts_json": section_alts_json,
            "groups_json": groups_json,
        }

        with open(out_path, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            w.writerow(row)

        return out_path
    
    def _write_ga_results_csv(self, out_path: Path, reports: list, ga_params: dict):
        """
        Write ALL GA candidate results to one CSV file.
        `reports` is the list returned by GAmod.GA(...).
        `ga_params` is a dict of the parameters used for the run (for reproducibility).
        """
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)

        # Define columns
        fieldnames = [
            "rank",
            "sequence",
            "structure",
            "fitness",
            "length",
            "melting_temp",
            "alt_structures",
            "placeholder",
            "placeholder_tm",
            "mirna_tm",
            "ga_params_json",
            "timestamp",
        ]

        ts = datetime.now().isoformat(timespec="seconds")
        ga_params_json = json.dumps(ga_params, ensure_ascii=False)

        with open(out_path, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            for rep in (reports or []):
                w.writerow({
                    "rank": rep.get("rank"),
                    "sequence": rep.get("sequence"),
                    "structure": rep.get("structure"),
                    "fitness": rep.get("fitness"),
                    "length": rep.get("length"),
                    "melting_temp": rep.get("melting_temp"),
                    "alt_structures": rep.get("alt_structures"),
                    "placeholder": rep.get("placeholder", ""),
                    "placeholder_tm": rep.get("placeholder_tm"),
                    "mirna_tm": rep.get("mirna_tm"),
                    "ga_params_json": ga_params_json,
                    "timestamp": ts,
                })




    # ----- Export -----
    def _export_all(self):
        if self.n == 0:
            messagebox.showwarning("Nothing to export", "Please render a structure first.")
            return

        desired = self.desired_structure

        # ----- Build seq_constraint with Section Alts only -----
        # Walk left-to-right; when a section alt starts, emit [ORIG|ALT1|ALT2] and skip its covered bases.
        sections = sorted(
            [{"start": max(0, s["start"]),
            "end":   min(self.n-1, s["end"]),
            "alts":  s["alts"]}
            for s in self.section_alts],
            key=lambda z: (z["start"], z["end"])
        )

        out_parts = []
        cursor = 0
        sec_idx = 0

        while cursor < self.n:
            took_section = False
            while sec_idx < len(sections) and sections[sec_idx]["end"] < cursor:
                sec_idx += 1
            if sec_idx < len(sections):
                sec = sections[sec_idx]
                if sec["start"] == cursor:
                    # ORIGINAL is whatever the left pane currently defines for those bases.
                    orig = "".join(self._label_for_allowed(self.allowed[k]) for k in range(sec["start"], sec["end"] + 1))
                    out_parts.append("[" + "|".join([orig] + sec["alts"]) + "]")
                    cursor = sec["end"] + 1
                    took_section = True
            if took_section:
                continue

            # No section here: emit per-base token from allowed set
            sset = self.allowed[cursor]
            if sset == set("ACGT"):
                out_parts.append("N")
            elif len(sset) == 1:
                out_parts.append(next(iter(sset)))
            else:
                out_parts.append("[" + "|".join(sorted(sset)) + "]")

            cursor += 1

        seq_str = "".join(out_parts)

        # ----- struct_constraints_base (digits) -----
        struct_digits = "".join(str(int(min(9, max(0, w)))) for w in self.weights)

        # ----- dynamic_penalty from groups (incl. terminals) -----
        dynamic_penalty = self._build_dynamic_penalty(struct_digits)

        # ----- Output -----
        self.out_desired.delete("1.0", tk.END); self.out_desired.insert("1.0", desired)
        self.out_seq.delete("1.0", tk.END); self.out_seq.insert("1.0", seq_str)
        self.out_struct.delete("1.0", tk.END); self.out_struct.insert("1.0", struct_digits)
        self.out_dyn.delete("1.0", tk.END); self.out_dyn.insert("1.0", dynamic_penalty)

        # ----- Save final export to CSV in Downloads -----
        try:
            out_path = self._write_export_csv(desired, seq_str, struct_digits, dynamic_penalty)
            self.status.set(f"Exported + saved CSV to: {out_path}")
        except Exception as e:
            # Don't block export if saving fails
            self.status.set(f"Exported (CSV save failed: {e})")


        self.status.set("Exported. Section alts are encoded as [ORIGINAL|ALT1|ALT2].")
        self.run_ga_btn.configure(state="normal")




if __name__ == "__main__":
    app = ConstraintUI()
    app.mainloop()
