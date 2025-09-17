import sys
import pandas as pd
import numpy as np
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
    QPushButton, QLabel, QLineEdit, QComboBox, QFileDialog, QTabWidget, QGroupBox,
    QCheckBox
)
from PyQt6.QtGui import QFont
from PyQt6.QtCore import QTimer, Qt
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import USB_1608FS_PLUS as daq


class USB1608FSPLUS(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('DAQ System - OM-USB-1608FS-PLUS')
        # Core DAQ Parameters
        self.sampling_rate = 50000
        self.buffer_size = 50000
        self.ai_range = ['BIP10VOLTS', 'BIP5VOLTS']
        self.daqs = ['OM-USB-1608FS-PLUS', 'NI']
        self.actual_sampling_rate = 0.0
        self.low_chan = 0
        self.high_chan = 0
        self.last_timestamp = 0.0
        # Data Buffers
        self.data_r = np.empty((0, 2))
        self.ylim = -5
        self.dist = 1
        self.ymax = 5
        # UI Logic Flags & Variables
        self._is_updating_channels = False
        self.selected_points = []
        self.markers = []
        self.annotations = []
        self.dragging_annotation = None
        self.drag_offset = (0, 0)
        # Plotting Objects
        self.signal_figure = None
        self.signal_canvas = None
        self.signal_axes = []
        self.active_plot_channels = []
        self.plot_map = {}
        self.plot_channel_config = {}
        # Main Timer
        self.update_timer = QTimer(self)
        self.update_timer.timeout.connect(self.read_and_update_plot)
        # Initialize UI
        self.setup_ui()

    def setup_ui(self):
        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        main_layout = QHBoxLayout(main_widget)
        # == Left Panel (Controls) ==
        left_panel = QWidget()
        left_panel.setFixedWidth(400)  # Increased width for the matrix
        left_panel_layout = QVBoxLayout(left_panel)
        left_panel_layout.setContentsMargins(10, 10, 10, 10)
        setup_group = QGroupBox("Setup")
        setup_group.setFont(QFont("Arial", 12, QFont.Weight.Bold))
        setup_grid = QGridLayout()
        self.board_num_en = QLineEdit(str(0))
        self.sampling_rate_en = QLineEdit(str(self.sampling_rate))
        self.ymin = QLineEdit(str(self.ylim))
        self.ymax = QLineEdit(str(self.ymax))
        self.dist = QLineEdit(str(self.dist))
        self.buffer_size_en = QLineEdit(str(self.buffer_size))
        self.daq_device_info = QComboBox()
        self.daq_device_info.addItems(self.daqs)
        self.channel_combo = QComboBox()
        self.window = QLineEdit(str(10))
        self.ai_range_en = QComboBox()
        self.ai_range_en.addItems(self.ai_range)
        setup_grid.addWidget(QLabel("Board Number:"), 0, 0)
        setup_grid.addWidget(self.board_num_en, 1, 0)
        setup_grid.addWidget(QLabel("Sampling Rate (S/s):"), 0, 1)
        setup_grid.addWidget(self.sampling_rate_en, 1, 1)
        setup_grid.addWidget(QLabel("Buffer Size:"), 2, 0)
        setup_grid.addWidget(self.buffer_size_en, 3, 0)
        setup_grid.addWidget(QLabel("DAQ Device:"), 4, 0)
        setup_grid.addWidget(self.daq_device_info, 5, 0)
        setup_grid.addWidget(QLabel("Window:"), 4, 1)
        setup_grid.addWidget(self.window, 5, 1)
        setup_grid.addWidget(QLabel("Analog Input Range:"), 2, 1)
        setup_grid.addWidget(self.ai_range_en, 3, 1)
        setup_group.setLayout(setup_grid)
        left_panel_layout.addWidget(setup_group)
        # Buttons GroupBox
        button_group = QGroupBox("Controls")
        button_group.setFont(QFont("Arial", 12, QFont.Weight.Bold))
        button_grid = QGridLayout()
        self.calibrate_btn = QPushButton("Calibration")
        self.calibrate_btn.setObjectName("calibrate_btn")
        self.start_btn = QPushButton("Start")
        self.start_btn.setObjectName("start_btn")
        self.reset_btn = QPushButton("Reset")
        self.reset_btn.setObjectName("reset_btn")
        self.stop_btn = QPushButton("Stop")
        self.stop_btn.setObjectName("stop_btn")
        self.save_btn = QPushButton("Save Data")
        self.save_btn.setObjectName("save_btn")
        self.stop_btn.setEnabled(False)
        self.reset_btn.setEnabled(False)
        button_grid.addWidget(self.calibrate_btn, 0, 0)
        button_grid.addWidget(self.start_btn, 0, 1)
        button_grid.addWidget(self.reset_btn, 1, 0)
        button_grid.addWidget(self.stop_btn, 1, 1)
        button_grid.addWidget(self.save_btn, 2, 0, 1, 2)  # Span across 2 columns
        button_group.setLayout(button_grid)
        left_panel_layout.addWidget(button_group)
        # --- TABS WIDGET ---
        self.control_tabs = QTabWidget()
        # ===================================================================
        # --- TAB 1: DAQ SETUP ---
        # ===================================================================
        self.daq_setup_tab = QWidget()
        daq_setup_layout = QVBoxLayout(self.daq_setup_tab)
        # --- Group 2: Channel Configuration Matrix ---
        channel_group = QGroupBox("Channel Configuration")
        channel_grid = QGridLayout(channel_group)
        headers = ["Active", "Channel", "Actual Rate (S/s)", "Buffer Size"]
        for col, text in enumerate(headers):
            label = QLabel(text);
            label.setFont(QFont("Arial", weight=QFont.Weight.Bold))
            channel_grid.addWidget(label, 0, col, Qt.AlignmentFlag.AlignCenter)
        self.channel_checkboxes = []
        self.channel_labels = []
        self.channel_rate_displays = []
        self.channel_buffer_displays = []
        for i in range(8):
            checkbox = QCheckBox()
            channel_label = QLabel(str(i))
            rate_display = QLineEdit();
            rate_display.setReadOnly(True);
            rate_display.setStyleSheet("background-color: #f0f0f0;")
            buffer_display = QLineEdit();
            buffer_display.setReadOnly(True);
            buffer_display.setStyleSheet("background-color: #f0f0f0;")
            self.channel_checkboxes.append(checkbox)
            self.channel_labels.append(channel_label)
            self.channel_rate_displays.append(rate_display)
            self.channel_buffer_displays.append(buffer_display)
            channel_grid.addWidget(checkbox, i + 1, 0, Qt.AlignmentFlag.AlignCenter)
            channel_grid.addWidget(channel_label, i + 1, 1, Qt.AlignmentFlag.AlignCenter)
            channel_grid.addWidget(rate_display, i + 1, 2)
            channel_grid.addWidget(buffer_display, i + 1, 3)

        daq_setup_layout.addWidget(channel_group)
        daq_setup_layout.addStretch()
        self.control_tabs.addTab(self.daq_setup_tab, "DAQ Setup")
        # ===================================================================
        # --- TAB 2: CALIBRATION ---
        # ===================================================================
        self.graph_tab = QWidget()
        graph_layout = QVBoxLayout(self.graph_tab)
        # --- Primer Grupo ---
        self.graph_group = QGroupBox("Graph Parameters")
        channel_graph = QGridLayout(self.graph_group)
        self.group_channels_checkbox = QCheckBox("Plotting Channels")
        channel_graph.addWidget(self.group_channels_checkbox, 0, 0)
        # --- Segundo Grupo ---
        self.scales_group = QGroupBox("Graph Y-Scale")
        scales_graph = QGridLayout(self.scales_group)
        # Añades widgets al segundo grupo
        scales_graph.addWidget(QLabel("Set Y min:"), 0, 0)
        scales_graph.addWidget(self.ymin, 1, 0)
        scales_graph.addWidget(QLabel("Set Y max:"), 0, 1)
        scales_graph.addWidget(self.ymax, 1, 1)
        # --- Terce Grupo ---
        # --- Segundo Grupo ---
        self.ponits_group = QGroupBox("Compute Points")
        ponits_graph = QGridLayout(self.ponits_group)
        # Añades widgets al segundo grupo
        ponits_graph.addWidget(QLabel("Set Dist:"), 0, 0)
        ponits_graph.addWidget(self.dist, 1, 0)
        # --------------------
        graph_layout.addWidget(self.graph_group)
        graph_layout.addWidget(self.scales_group)
        graph_layout.addWidget(self.ponits_group)
        graph_layout.addStretch()
        self.control_tabs.addTab(self.graph_tab, "Graph Config")
        left_panel_layout.addWidget(self.control_tabs)
        left_panel_layout.addWidget(button_group)
        # == Right Panel (Plots) ==
        self.tabs = QTabWidget()
        self.signal_tab = QWidget()
        self.signal_tab_layout = QVBoxLayout(self.signal_tab)
        self.tabs.addTab(self.signal_tab, "Signal Analysis")
        # == Final Assembly & Connections ==
        main_layout.addWidget(left_panel)
        main_layout.addWidget(self.tabs)
        self.start_btn.clicked.connect(self.start_acquisition)
        self.stop_btn.clicked.connect(self.stop_acquisition)
        self.reset_btn.clicked.connect(self.reset_data)
        self.save_btn.clicked.connect(self.save_data)
        self.sampling_rate_en.textChanged.connect(self._calculate_and_update_actual_rate)
        self.buffer_size_en.textChanged.connect(self._update_buffer_displays)
        for checkbox in self.channel_checkboxes:
            checkbox.stateChanged.connect(self._update_cascading_channels)
        # Set the initial UI state
        self.channel_checkboxes[0].setChecked(True)

    def _set_row_visibility(self, row_index, is_visible):
        self.channel_checkboxes[row_index].setVisible(is_visible)
        self.channel_labels[row_index].setVisible(is_visible)
        self.channel_rate_displays[row_index].setVisible(is_visible)
        self.channel_buffer_displays[row_index].setVisible(is_visible)

    def _update_cascading_channels(self):
        if self._is_updating_channels: return
        self._is_updating_channels = True
        first_unchecked_index = -1
        for i, cb in enumerate(self.channel_checkboxes):
            if not cb.isChecked():
                first_unchecked_index = i
                break
        for i in range(8):
            if first_unchecked_index == -1 or i <= first_unchecked_index:
                self._set_row_visibility(i, True)
                if first_unchecked_index == -1 or i < first_unchecked_index:
                    self.channel_checkboxes[i].setChecked(True)
                else:
                    self.channel_checkboxes[i].setChecked(False)
            else:
                self._set_row_visibility(i, False)
                self.channel_checkboxes[i].setChecked(False)
        self._is_updating_channels = False
        self._update_plots_layout()
        self._calculate_and_update_actual_rate()

    def _calculate_and_update_actual_rate(self):
        MAX_AGGREGATE_RATE = 400000.0
        try:
            desired_rate = float(self.sampling_rate_en.text())
        except ValueError:
            desired_rate = 0
        active_flags = [cb.isChecked() for cb in self.channel_checkboxes]
        num_active = sum(active_flags)
        actual_rate = 0
        color = "#f0f0f0"
        if num_active > 0:
            if desired_rate * num_active > MAX_AGGREGATE_RATE:
                actual_rate = MAX_AGGREGATE_RATE / num_active
                color = "#fff3cd"  # Yellow
            else:
                actual_rate = desired_rate
                color = "#d4edda"  # Green
        rate_text = f"{actual_rate:.2f}"
        for i, field in enumerate(self.channel_rate_displays):
            if active_flags[i]:
                field.setText(rate_text)
                field.setStyleSheet(f"background-color: {color};")
            else:
                field.setText("")
                field.setStyleSheet("background-color: #f0f0f0;")
        self._update_buffer_displays()

    def _update_buffer_displays(self):
        buffer_text = self.buffer_size_en.text()
        active_flags = [cb.isChecked() for cb in self.channel_checkboxes]
        if hasattr(self, 'channel_buffer_displays'):
            for i, field in enumerate(self.channel_buffer_displays):
                if active_flags[i]:
                    field.setText(buffer_text)
                else:
                    field.setText("")

    def _update_plots_layout(self):
        self._clean_plot_widgets()
        self.active_plot_channels = [i for i, cb in enumerate(self.channel_checkboxes) if cb.isChecked()]
        self.plot_map = {}
        self.plot_channel_config = {}
        active_groups = []

        if self.group_channels_checkbox.isChecked():
            if self.active_plot_channels:
                active_groups = ['all_in_one']
                self.plot_channel_config['all_in_one'] = self.active_plot_channels
        else:
            default_channel_groups = {0: [0, 4], 1: [1, 5], 2: [2, 6], 3: [3, 7]}
            active_groups = [
                group for group, channels in default_channel_groups.items()
                if any(ch in self.active_plot_channels for ch in channels)
            ]
            self.plot_channel_config = default_channel_groups
        self.signal_figure = Figure(constrained_layout=True)
        self.signal_canvas = FigureCanvas(self.signal_figure)
        self.toolbar = NavigationToolbar(self.signal_canvas, self)
        self.signal_tab_layout.addWidget(self.toolbar)
        self.signal_tab_layout.addWidget(self.signal_canvas)
        self.signal_axes = []
        if not active_groups:
            self.signal_canvas.draw()
            return
        for i, group_key in enumerate(active_groups):
            ax = self.signal_figure.add_subplot(len(active_groups), 1, i + 1)
            self.signal_axes.append(ax)
            self.plot_map[group_key] = ax
        if len(active_groups) > 1:
            self.signal_figure.subplots_adjust(hspace=0.6)
        self.signal_canvas.mpl_connect('button_press_event', self.on_click)
        self.signal_canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.signal_canvas.mpl_connect('button_release_event', self.on_release)

    def _setup_new_acquisition(self):
        print("Performing full setup for new acquisition...")
        self._update_plots_layout()
        selected_channels = [i for i, cb in enumerate(self.channel_checkboxes) if cb.isChecked()]
        if not selected_channels:
            raise ValueError("No channels selected.")
        self.low_chan = min(selected_channels)
        self.high_chan = max(selected_channels)
        num_scanned = self.high_chan - self.low_chan + 1
        self.data_r = np.empty((0, 1 + num_scanned))

    def _initialize_data_buffer(self):
        selected_channels = [i for i, cb in enumerate(self.channel_checkboxes) if cb.isChecked()]
        if selected_channels:
            self.low_chan = min(selected_channels)
            self.high_chan = max(selected_channels)
            num_scanned = self.high_chan - self.low_chan + 1
            self.data_r = np.empty((0, 1 + num_scanned))
        else:
            self.data_r = np.empty((0, 2))

        print(f"Data buffer initialized with {self.data_r.shape[1]} columns.")

    def start_acquisition(self):
        selected_channels = [i for i, cb in enumerate(self.channel_checkboxes) if cb.isChecked()]
        if not selected_channels:
            print("Error: Select at least one channel.")
            self.set_controls_enabled(True)
            return
        num_new_chans = max(selected_channels) - min(selected_channels) + 1
        config_changed = (self.data_r.ndim < 2) or (num_new_chans != self.data_r.shape[1] - 1)
        if self.signal_canvas is None or config_changed:
            self.stop_acquisition()
            self._clean_plot_widgets()
            self._update_plots_layout()
            self._initialize_data_buffer()
            daq.reset_time_counter()
        self.set_controls_enabled(False)
        try:
            board_num = int(self.board_num_en.text())
            rate = int(self.sampling_rate_en.text())
            samples_per_chan = int(self.buffer_size_en.text())
            range_text = self.ai_range_en.currentText()
            voltage_range = getattr(daq, range_text)
            num_scanned = self.high_chan - self.low_chan + 1
            daq.start_stream(board_num, self.low_chan, self.high_chan,
                             samples_per_chan, rate, voltage_range)
            self.board_num = board_num
            self.actual_sampling_rate = min(400000.0 / num_scanned, float(rate))
            self.update_timer.start(100)

        except Exception as e:
            print(f"Error starting acquisition: {e}")
            self._clean_plot_widgets()
            self.set_controls_enabled(True)

    def stop_acquisition(self):
        self.update_timer.stop()
        if hasattr(self, 'data_r') and self.data_r.shape[0] > 0:
            self.last_timestamp = self.data_r[-1, 0]
        if hasattr(self, 'board_num'):
            daq.stop_stream(self.board_num)
        self.set_controls_enabled(True)

    def _clean_plot_widgets(self):
        for marker in self.markers:
            try:
                marker.remove()
            except:
                pass
        for annotation in self.annotations:
            try:
                annotation.remove()
            except:
                pass
        self.markers.clear()
        self.annotations.clear()
        self.selected_points.clear()

        if hasattr(self, 'toolbar') and self.toolbar:
            self.signal_tab_layout.removeWidget(self.toolbar)
            self.toolbar.deleteLater()
            self.toolbar = None

        if hasattr(self, 'signal_canvas') and self.signal_canvas:
            self.signal_tab_layout.removeWidget(self.signal_canvas)
            self.signal_canvas.deleteLater()
            self.signal_canvas = None

    def read_and_update_plot(self):
        try:
            new_data = daq.read_block(self.board_num)
            if new_data.shape[0] > 0:
                self.data_r = np.vstack([self.data_r, new_data])
            if not hasattr(self, 'plot_map') or self.data_r.shape[0] == 0:
                return
            window_seconds = int(self.window.text())
            samples_to_show = int(self.actual_sampling_rate * window_seconds)
            plot_data = self.data_r[-samples_to_show:, :]
            time_array = plot_data[:, 0]
            plot_colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
                           "#9467bd", "#8c564b", "#e377c2", "#7f7f7f"]
            for ax in self.signal_axes:
                ax.clear()
            for group_key, ax in self.plot_map.items():
                channels_to_plot = self.plot_channel_config.get(group_key, [])
                active_in_group = []
                for ch in channels_to_plot:
                    if ch in self.active_plot_channels:
                        active_in_group.append(ch)
                        data_col = 1 + (ch - self.low_chan)
                        if data_col < plot_data.shape[1]:
                            voltage = plot_data[:, data_col]
                            ax.plot(time_array, voltage,
                                    label=f'Ch{ch}',
                                    linewidth=0.8,
                                    color=plot_colors[ch % len(plot_colors)])

                if active_in_group:
                    ax.legend(loc='upper right')
                    ax.grid(False)
                    ax.set_ylim(float(self.ymin.text()), float(self.ymax.text()))

            if self.signal_axes:
                for ax in self.signal_axes:
                    ax.tick_params(labelbottom=True)
                self.signal_axes[-1].set_xlabel('Time (s)')

            self.signal_canvas.draw_idle()

        except Exception as e:
            print(f"Error during plot update: {e}")
            self.stop_acquisition()

    def set_controls_enabled(self, is_enabled):
        for widget in [self.sampling_rate_en, self.buffer_size_en, self.daq_device_info, self.ai_range_en,
                       self.start_btn, self.reset_btn]:
            widget.setEnabled(is_enabled)
        for cb in self.channel_checkboxes:
            cb.setEnabled(is_enabled)
        self.stop_btn.setEnabled(not is_enabled)

    def reset_data(self):
        self.stop_acquisition()
        self._clean_plot_widgets()
        self.data_r = np.empty((0, 2))

        daq.reset_time_counter()

        self.set_controls_enabled(True)
        self.reset_btn.setEnabled(False)

    def save_data(self):
        if self.data_r.shape[0] == 0:
            print("No data to save.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Save Data", "", "CSV Files (*.csv)")
        if path:
            try:
                num_chans = self.data_r.shape[1] - 1
                headers = ['Time (s)'] + [f'Channel {self.low_chan + i}' for i in range(num_chans)]
                df = pd.DataFrame(self.data_r, columns=headers)
                df.to_csv(path, index=False, float_format='%.6f')
                print(f"Data saved to {path}")
            except Exception as e:
                print(f"Failed to save data: {e}")

    def on_click(self, event):
        if event.inaxes is None or self.signal_canvas.toolbar.mode != '':
            return
        ax = event.inaxes
        if event.button == 3:
            self.clear_selections(ax)
            self.signal_canvas.draw_idle()
            return
        if event.button != 1:
            return
        for ann in self.annotations:
            if ann.axes == ax:
                contains, _ = ann.contains(event)
                if contains:
                    self.dragging_annotation = ann
                    self.drag_offset = (ann.xyann[0] - event.x, ann.xyann[1] - event.y)
                    return
        channels_in_ax = []
        for group_key, group_ax in self.plot_map.items():
            if group_ax == ax:
                config = self.plot_channel_config.get(group_key, [])
                channels_in_ax = [ch for ch in config if ch in self.active_plot_channels]
                break
        if not channels_in_ax or self.data_r.shape[0] == 0:
            return
        min_vertical_dist = float('inf')
        best_match = None
        window_seconds = int(self.window.text())
        samples_to_show = int(self.actual_sampling_rate * window_seconds)
        plot_data = self.data_r[-samples_to_show:, :]
        time_array = plot_data[:, 0]
        for ch in channels_in_ax:
            data_col = 1 + (ch - self.low_chan)
            if data_col >= plot_data.shape[1]:
                continue
            voltage_data = plot_data[:, data_col]
            idx = (np.abs(time_array - event.xdata)).argmin()
            x_real, y_real = time_array[idx], voltage_data[idx]
            vertical_dist = abs(y_real - event.ydata)
            if vertical_dist < min_vertical_dist:
                min_vertical_dist = vertical_dist
                best_match = {'x': x_real, 'y': y_real, 'channel': ch}
        if best_match:
            punto = ax.plot(best_match['x'], best_match['y'], 'ro', markersize=5)[0]
            self.markers.append(punto)
            point_info = {
                'x': best_match['x'],
                'y': best_match['y'],
                'ax': ax,
                'channel': best_match['channel']
            }
            self.selected_points.append(point_info)
        if len(self.selected_points) == 2:
            p1, p2 = self.selected_points
            if p1['ax'] == p2['ax']:
                linea = p1['ax'].plot([p1['x'], p2['x']], [p1['y'], p2['y']], 'k--', linewidth=1)[0]
                self.markers.append(linea)
            dx = abs(p2['x'] - p1['x'])
            dy = abs(p2['y'] - p1['y'])
            texto = (f"P1 (Ch{p1['channel']}): ({p1['x']:.6f}s, {p1['y']:.6f}V)\n"
                     f"P2 (Ch{p2['channel']}): ({p2['x']:.6f}s, {p2['y']:.6f}V)\n"
                     f"Δt={dx:.6f}s, ΔV={dy:.6f}V\n"
                     f"Velocity = {float(self.dist.text()) / dx:.6f} m/s")
            anotacion = p1['ax'].annotate(
                texto, xy=(0.02, 0.98), xycoords='axes fraction',
                textcoords="offset pixels", xytext=(0, -5),
                fontsize=9, va='top', ha='left',
                bbox=dict(boxstyle="round,pad=0.4", fc="lightyellow", alpha=0.8)
            )
            self.annotations.append(anotacion)
            self.selected_points.clear()
        self.signal_canvas.draw_idle()
    def on_motion(self, event):
        if self.dragging_annotation is None or event.inaxes != self.dragging_annotation.axes:
            return
        new_xy = (event.x + self.drag_offset[0], event.y + self.drag_offset[1])
        self.dragging_annotation.xyann = new_xy
        self.signal_canvas.draw_idle()
    def on_release(self, event):
        self.dragging_annotation = None
    def clear_selections(self, ax):
        markers_a_borrar = [m for m in self.markers if m.axes == ax]
        for m in markers_a_borrar:
            m.remove()
            self.markers.remove(m)
        annotations_a_borrar = [a for a in self.annotations if a.axes == ax]
        for a in annotations_a_borrar:
            a.remove()
            self.annotations.remove(a)
        self.selected_points = [p for p in self.selected_points if p['ax'] != ax]
def main():
    app = QApplication(sys.argv)
    style_sheet = """
        /* Estilo base para todos los botones */
        QPushButton {
            color: white;
            padding: 6px;
            border-radius: 5px;
            font-weight: bold;
            border: 1px solid rgba(0, 0, 0, 0.2);
        }

        /* 🟢 Botones HABILITADOS (verdes) */
        QPushButton:enabled {
            background-color: #2ecc71; /* verde */
        }
        QPushButton:enabled:hover {
            background-color: #27ae60; /* verde más oscuro */
        }

        /* 🔴 Botones INHABILITADOS (rojos) */
        QPushButton:disabled {
            background-color: #e74c3c; /* rojo */
            color: #ecf0f1;
        }

        /* 🔵 Botón de calibración específico (azul) - SOBREESCRIBE los demás */
        QPushButton#calibrate_btn {
            background-color: #3498db; /* azul */
        }
        QPushButton#calibrate_btn:hover {
            background-color: #2980b9; /* azul más oscuro */
        }

        /* Estilos para otros widgets */
        QLineEdit, QComboBox {
            padding: 3px;
            border: 1px solid gray;
            border-radius: 4px;
        }
    """
    app.setStyleSheet(style_sheet)
    window = USB1608FSPLUS()
    window.showMaximized()
    sys.exit(app.exec())


if __name__ == '__main__':
    main()
