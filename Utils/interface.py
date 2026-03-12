"""
    This file gather all functions and classes related to the user interface design.

    @author: Maxime Bernard

"""

from PyQt5.QtWidgets import (QErrorMessage,QSplitter,QTabWidget,QMainWindow,
                             QMdiArea,QDesktopWidget,QMessageBox,QMdiSubWindow,
                             QProgressBar,QLabel,QHBoxLayout,QPushButton,
                             QVBoxLayout,QPlainTextEdit,QWidget, QFileDialog,
                             QInputDialog,QLineEdit,QToolBar,QAction,QDockWidget,
                             QCheckBox,QComboBox,QGroupBox,QGridLayout,QTableWidget,
                             QTableWidgetItem,QApplication,qApp,QTabBar,QStylePainter,
                             QStyleOptionTab,QStyle,QFrame)
from PyQt5.QtGui import (QIcon,QKeySequence,QFont)
import configs as conf
import os
from PyQt5.QtCore import (Qt)
from PyQt5 import QtCore
import csv
import io
import matplotlib.patches as patches
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import numpy as np

#----------------------------------------------------------
def createTable(d, c):
    """ To create table as input parameter."""
    
    table = QTableWidget()
    table.setRowCount(d)
    table.setColumnCount(c)
    return table

#----------------------------------------------------------
class emptyText(QWidget):
    """ 
     This class stores the name of the Pecube project.
     
    """
    def __init__(self):
        super().__init__()
        self.ProjectName = ""

#----------------------------------------------------------
def progressbar(output,progress_re,p):
    """ Define a progress bar when running a Pecube model"""
    m = progress_re.search(output)
    if m:
        if p==1:
            pc_complete = 100-int(m.group(3))
        elif p==2:
            pc_complete = int(m.group(1))/int(m.group(2))*100
        return int(pc_complete)
    


##############################################################################
############################## Classes #######################################
##############################################################################


# Classes related to the graphical interface 
#############################################################
class QCustomSplitter(QSplitter):
    """ 
    A custom class for the splitters. A splitter if a moveable delimitation
    between two space in a window 
     
     """
    def __init__(self, parent=None):
        super(QCustomSplitter,self).__init__(parent)
        self.setStyleSheet("QSplitterHandle:hover {}  QSplitter::handle:hover {background-color:rgba(150, 150, 150, 100%);}")

    
#############################################################   
class QCustomTabWidget(QTabWidget):
    """ 
    A custom class for tabs. Enable the user to close the tab by calling
    the function closeTab 
    
    """
    def __init__(self, parent=None):
        super(QCustomTabWidget, self).__init__(parent)
        self.setTabsClosable(True)
        self.tabCloseRequested.connect(self.closeTab)
        
    #----------------------------------------------------------
    def closeTab(self, currentIndex):
        self.setTabVisible(currentIndex,False)   
        
        
#############################################################   
class run_Window(QMainWindow):   
    """ 
    This class open a window to embed sub-windows for each Pecube project run
    
    See also: MainWindow
    """
    
    def __init__(self, parent, height=800, width=1100,text='Running models...',fixed=False):
        super().__init__()
        self.mdi = QMdiArea()
        self.setCentralWidget(self.mdi)
        self.parent = parent
        self.closeAll = 0 # Signal to close all windows
        
        self.setWindowTitle(text)
        self.setWindowIcon(QIcon(parent.IconPecube))
        self.setGeometry(800, 400, width, height)
        if fixed:
            self.setMaximumWidth(width)
            self.setMaximumHeight(height)
        #force the window to open at the center of the screen     
        qtRectangle = self.frameGeometry()
        centerPoint = QDesktopWidget().availableGeometry().center()
        qtRectangle.moveCenter(centerPoint)
        self.move(qtRectangle.topLeft())
        
    def closeEvent(self, event):
        #to ask the user to confirm the closure of all running windows
        reply = QMessageBox.question(self, "Window Close", "Are you sure you want to close this window?"+
                                     "\n All running models will be stopped...",
                                      QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if reply == QMessageBox.Yes:
            self.closeAll = 1
            for i in self.mdi.subWindowList():
                if i.parent.p: # if a process is running
                    i.parent.p.kill() #kill process
                    i.close()
            self.parent.runOpen = 0
            #Clear from Window open list
            try:
                for i in range(len(conf.WindowsOpen)):
                    if conf.WindowsOpen[i] == self:
                        conf.WindowsOpen.pop(i)
                        i = i-1
            except IndexError:
                pass
            event.accept()
        else:
            event.ignore()
            
      
#############################################################   
class showMessage(QMdiSubWindow):
    """ This class is used to show the Pecube run outputs in a sub-window.
    
    console is the window where the output messages for the run are written
    mainWin is the main window object
    
    See also: run_Window
    """
    
    def __init__(self,parent, console, mainWin, ProjectName):
        super().__init__()

        self.parent = parent
        self.mainWin = mainWin
        self.progress = QProgressBar(self)
        self.progress.setRange(0, 100)
        self.progress.setValue(0)
        self.Label = QLabel('Pecube is running...')
        self.Label.setFont(conf.fontBold12)
        self.Label.setAlignment(Qt.AlignCenter)
        hbox = QHBoxLayout()
        self.OkButton = QPushButton('Ok')
        self.OkButton.setEnabled(False)
        self.CancelButton = QPushButton('Cancel')
        self.OkButton.clicked.connect(lambda: self.close())
        hbox.addStretch(1)
        hbox.addWidget(self.OkButton)
        hbox.addWidget(self.CancelButton)
        hbox.addStretch(1)
        
        IconWindow = os.path.join(conf.IconPath,"Pecube_icon.ico")
        self.setWindowIcon(QIcon(IconWindow))
        self.setWindowTitle(ProjectName)
        self.setGeometry(800, 400, 300, 200)
        #force the window to open at the center of the screen     
        qtRectangle = self.frameGeometry()
        centerPoint = QDesktopWidget().availableGeometry().center()
        qtRectangle.moveCenter(centerPoint)
        self.move(qtRectangle.topLeft())
        
        l = QVBoxLayout()
        l.addStretch(1)
        l.addWidget(self.Label)
        l.addWidget(self.progress)
        #if Console, show also text output from Pecube
        if int(console) == 1:
            self.resize(500,400)
            self.text = QPlainTextEdit()
            self.text.setReadOnly(True)
            l.addWidget(self.text,5)
        l.addStretch(1)
        l.addLayout(hbox)
        
        w = QWidget()
        w.setLayout(l)

        self.setWidget(w)
        
        # add sub-window to object "run_Window"
        mainWin.runWindow.mdi.addSubWindow(self)
        self.show()
        
    #------------------------------------------------------------------
    def closeEvent(self, event):
        #to ask the user to confirm the closure of the runing window
        if self.parent.p and not self.mainWin.runWindow.closeAll:
            reply = QMessageBox.question(self, "Window Close", "Are you sure you want to close this window?"+
                                     "\n The current running model will be stopped...",
                                      QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        
            if reply == QMessageBox.Yes:
                try:
                    self.parent.p.kill() #kill process
                except AttributeError:
                    print("process already canceled.")
                event.accept()
                self.parent.aborded = 1
            else:
                event.ignore()
                
        elif self.mainWin.runWindow.closeAll:
            event.accept()
            self.parent.aborded = 1
        self.parent.RunWindow = 0


#############################################################
class HorizontalTabWidget(QTabBar):
    """ To make the text horizontal in the tabs of input parameters. """
    def paintEvent(self, event):
        painter = QStylePainter(self)
        option = QStyleOptionTab()
        for index in range(self.count()):
            self.initStyleOption(option, index)
            painter.drawControl(QStyle.CE_TabBarTabShape, option)
            painter.drawText(self.tabRect(index),
                              QtCore.Qt.AlignCenter | QtCore.Qt.TextDontClip,
                              self.tabText(index))
            
    #----------------------------------------------------------
    def tabSizeHint(self, index):
        size = QTabBar.tabSizeHint(self, index)
        if size.width() < size.height():
            size.transpose()
        return size
    
    
#############################################################
class WinTable(QTableWidget):
    """
    Table to feed with data.
    
    """
    
    def __init__(self, data, *args):
        QTableWidget.__init__(self, *args)
        self.data = data
        self.setData()

    #----------------------------------------------------------
    def setData(self):
        HeaderLabel = []
        for n, key in enumerate(self.data.keys()):
            HeaderLabel.append(key)
            for m, item in enumerate(self.data[key]):
                newItem = QTableWidgetItem(item)
                newItem.setFont(conf.font6)
                self.setItem(m, n, newItem)
        self.setHorizontalHeaderLabels(HeaderLabel)

    #----------------------------------------------------------
    def keyPressEvent(self, event):
        clipboard = QApplication.clipboard()
        
        if event.matches(QKeySequence.Copy) and (event.modifiers() & Qt.ControlModifier):
            # self.copied_cells = sorted(self.selectedIndexes())
            self.copied_cells = self.copySelection()
            try:
                clipboard.setText()
            except TypeError:
                pass
        elif event.matches(QKeySequence.Paste) and (event.modifiers() & Qt.ControlModifier):
            self.pasteItem()
            # r = self.currentRow() - self.copied_cells[0].row()
            # c = self.currentColumn() - self.copied_cells[0].column()
            # for cell in self.copied_cells:
            #     self.setItem(cell.row() + r, cell.column() + c, QTableWidgetItem(cell.data()))

    #----------------------------------------------------------
    def copySelection(self):
        #to be able to copy data in the tables
        selection = self.selectedIndexes()
        
        if selection:
            rows = sorted(index.row() for index in selection)
            columns = sorted(index.column() for index in selection)
            rowcount = rows[-1] - rows[0] + 1
            colcount = columns[-1] - columns[0] + 1
            table = [[''] * colcount for _ in range(rowcount)]
            for index in selection:
                row = index.row() - rows[0]
                column = index.column() - columns[0]
                table[row][column] = index.data()
            stream = io.StringIO()
            csv.writer(stream, delimiter='\t').writerows(table)
            qApp.clipboard().setText(stream.getvalue())
        return

    #----------------------------------------------------------
    def pasteItem(self):
        #to be able to paste data in the tables
        selection = self.selectedIndexes()
        if selection:
            model = self.model()
            buffer = qApp.clipboard().text()
            rows = sorted(index.row() for index in selection)
            columns = sorted(index.column() for index in selection)
            reader = csv.reader(io.StringIO(buffer), delimiter='\t')
            
            if len(rows) == 1 and len(columns) == 1:
                for i, line in enumerate(reader):
                    for j, cell in enumerate(line):
                        model.setData(model.index(
                            rows[0]+i, columns[0]+j), cell)
            else:
                arr = [[cell for cell in row] for row in reader]
                for index in selection:
                    row = index.row() - rows[0]
                    column = index.column() - columns[0]
                    model.setData(model.index(
                        index.row(), index.column()), arr[row][column])
        return
    
    
#############################################################
class WinExtraParameters(QMainWindow):
    """
      This class enables to open and defines a new window to set extra Pecube parameters.
      when calling this class you can set the width, height and the text of the
      window.
      
    """
    def __init__(self, width=400, height=600,text='Pecube extra parameters',fixed=False):
        super().__init__()
        #FolderPath = QtCore.QCoreApplication.applicationDirPath()
        IconStep4 = os.path.join(conf.IconPath,"Pecube_icon.ico")
        self.setWindowTitle(text)
        self.setWindowIcon(QIcon(IconStep4))
        self.setGeometry(400, 200, width,height)
        if fixed:
            self.setMaximumWidth(width)
            self.setMaximumHeight(height)
        
        #force the window to open at the center of the screen     
        qtRectangle = self.frameGeometry()
        centerPoint = QDesktopWidget().availableGeometry().center()
        qtRectangle.moveCenter(centerPoint)
        self.move(qtRectangle.topLeft())
        
        self.layout = QGridLayout()
        hbox = QHBoxLayout()
        self.OkButton = QPushButton('Ok')
        self.CancelButton = QPushButton('Cancel')
        # self.OkButton.clicked.connect(lambda: self.close())
        self.CancelButton.clicked.connect(lambda: self.close())
        hbox.addStretch(1)
        hbox.addWidget(self.OkButton)
        hbox.addWidget(self.CancelButton)
        hbox.addStretch(1)
        self.layout.addLayout(hbox, 3, 0, 1, 3)


#############################################################
class oldProjectName(QWidget):
    """
    
    Search for the old Pecube project name when opening an old input file
    The project Name is read from the path folder of the input file
    which is always *\ProjectName\Input\Pecube.in
    Because the Project Name has to be 5 letters long, it corresponds to
    the position in the folder path, -21:-16
    
    """
    def __init__(self,parent):
        super().__init__()
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        dialog = QFileDialog()
        dialog.setDirectory(os.path.join(parent.PecubePath))
        self.fileName, _ = dialog.getOpenFileName(
            self, "QFileDialog.getOpenFileName()", "", "Pecube Files (*.in)", options=options)
        self.ProjectName = self.fileName[-21:-16]
        

#########################################################################          
class DraggablePoint:
    def __init__(self, parent, enable, x=[0,1], y=[20,20], size=1, color ='r', ax=None):
        """Creates a draggable Point on a matplotlib canvas"""
        # The FigureCanvas
        self.parent = parent
        self.point_alpha_default = 0.8
        self.x = x
        self.n = 2
        self.size = size
        self.mousepress = None
        self.currently_dragging = False
        self.current_artist = None
        self.offset = [0,0]
        self.listLabelPoints = []
        self.color = color
        self.ax = ax
        self.enable = enable
        
        # Add points
        # Grain surface (=0)
        for i in range(len(self.x)):
            point_object = patches.Ellipse((x[i], y[i]), 0.015, size, fc=self.color, edgecolor=self.color,
                   alpha=self.point_alpha_default, transform=ax.transData, label="point"+str(i))
            point_object.set_picker(5)
            ax.add_patch(point_object)
            self.listLabelPoints.append('point'+str(i))
        
        xdata = []
        ydata = []
        for p in ax.patches:
            cx, cy = p.center
            xdata.append(cx)
            ydata.append(cy)
        self.line_object = ax.plot(xdata, ydata, alpha=0.5, c=self.color, lw=2, picker=True)
        self.line_object[0].set_pickradius(5)
        
        self.parent.draw()
        
        self.connect()
    
    #------------------------------------------------
    def connect(self):  
        self.parent.fig.canvas.mpl_connect('button_press_event', self.on_press)
        self.parent.fig.canvas.mpl_connect('button_release_event', self.on_release)
        self.parent.fig.canvas.mpl_connect('pick_event',self.on_pick)
        self.parent.fig.canvas.mpl_connect('motion_notify_event', self.on_motion)  
      
    #------------------------------------------------
    def on_press(self,event):
        self.currently_dragging = True
        if event.button == 3:
            self.mousepress = "right"
        elif event.button == 1:
            self.mousepress = "left"
    
    #------------------------------------------------
    def on_release(self,event):
        self.current_artist = None
        self.currently_dragging = False
    
    #------------------------------------------------
    def on_pick(self,event):
        
        #First check which data is picked (U or Th)
        if self.enable:
            if self.current_artist is None:
                self.current_artist = event.artist
                #print("pick ", current_artist)
                if isinstance(event.artist, patches.Ellipse):
                    if event.mouseevent.dblclick:
                        if self.mousepress == "right":
                            # print("double click right")
                            if len(self.ax.patches) > 2:
                                #print("\ndelete", event.artist.get_label())
                                event.artist.remove()
                                xdata = list(self.line_object[0].get_xdata())
                                ydata = list(self.line_object[0].get_ydata())
                                for i in range(0,len(xdata)):
                                    if event.artist.get_label() == self.listLabelPoints[i]:
                                        xdata.pop(i) 
                                        ydata.pop(i) 
                                        self.listLabelPoints.pop(i)
                                        break
                                #print('--->', listLabelPoints)
                                self.line_object[0].set_data(xdata, ydata)
                                self.parent.draw()
                    else:
                        x0, y0 = self.current_artist.center
                        x1, y1 = event.mouseevent.xdata, event.mouseevent.ydata
                        self.offset = [(x0 - x1), (y0 - y1)]
                        
                elif isinstance(event.artist, Line2D):
                    if event.mouseevent.dblclick:
                        if self.mousepress == "left":
                            # print("double click left")
                            self.n = self.n+1
                            x, y = event.mouseevent.xdata, event.mouseevent.ydata
                            newPointLabel = "point"+str(self.n)
                            point_object = patches.Ellipse((x, y), 0.015, self.size, fc=self.color, edgecolor=self.color,
                                   alpha=self.point_alpha_default, transform=self.ax.transData, label=newPointLabel)
                            point_object.set_picker(5)
                            xdata = list(self.line_object[0].get_xdata())
                            ydata = list(self.line_object[0].get_ydata())
                            #print('\ninit', listLabelPoints)
                            pointInserted = False
                            for i in range(0,len(xdata)-1):
                                #print("--> testing inclusion %s in [%s-%s]" 
                                #        %(newPointLabel, listLabelPoints[i], listLabelPoints[i+1]))
                                toly = 1;
                                # print('----->', min(xdata[i],xdata[i+1]), '<', x, '<', max(xdata[i],xdata[i+1]))
                                # print('----->', min(ydata[i],ydata[i+1])-toly, '<', y, '<', max(ydata[i],ydata[i+1])+toly)
                                if x > min(xdata[i],xdata[i+1]) and x < max(xdata[i],xdata[i+1]) and \
                                   y > min(ydata[i],ydata[i+1])-toly and y < max(ydata[i],ydata[i+1]) +toly:
                                    xdata.insert(i+1, x)
                                    ydata.insert(i+1, y)
                                    self.listLabelPoints.insert(i+1, newPointLabel)
                                    pointInserted = True
                                    self.ax.add_patch(point_object)
                                    break
                            if not pointInserted:
                                print("Error: point not inserted")
                                return
                            self.line_object[0].set_data(xdata, ydata)
                            #print('final', listLabelPoints)
                            self.parent.draw()    
                    else:
                        xdata = event.artist.get_xdata()
                        ydata = event.artist.get_ydata()
                        x1, y1 = event.mouseevent.xdata, event.mouseevent.ydata
                        self.offset = xdata[0] - x1, ydata[0] - y1
    
    #------------------------------------------------
    def on_motion(self,event):
        if not self.enable:
            return
        if not self.currently_dragging:
            return
        if self.current_artist == None:
            return
        if event.xdata == None:
            return
        dx,dy = self.offset
        if isinstance(self.current_artist, patches.Ellipse):
            if self.current_artist.get_label() == self.listLabelPoints[0]:#don't move edge point through x
                cx = 0.0
                cy = event.ydata + dy
            elif self.current_artist.get_label() == self.listLabelPoints[-1]:
                cx = 1.0
                cy = event.ydata + dy
            else:
                cx, cy =  event.xdata + dx, event.ydata + dy
            self.current_artist.center = cx, cy
            #print("moving", current_artist.get_label())
            xdata = list(self.line_object[0].get_xdata())
            ydata = list(self.line_object[0].get_ydata())
            for i in range(0,len(xdata)): 
                if self.listLabelPoints[i] == self.current_artist.get_label():
                    xdata[i] = cx
                    ydata[i] = cy
                    break
            self.line_object[0].set_data(xdata, ydata)
        elif isinstance(self.current_artist, Line2D):
            xdata = list(self.line_object[0].get_xdata())
            ydata = list(self.line_object[0].get_ydata())
            xdata0 = xdata[0]
            ydata0 = ydata[0]
            for i in range(0,len(xdata)):
                    xdata[i] = xdata[i] #Keep all points static on x axis
                    ydata[i] = event.ydata + dy + ydata[i] - ydata0 
            self.line_object[0].set_data(xdata, ydata)
            for p in self.ax.patches:
                pointLabel = p.get_label()
                i = self.listLabelPoints.index(pointLabel) 
                p.center = xdata[i], ydata[i]
        self.parent.draw()
  
    
#########################################################################      
class InteractivePlot(FigureCanvas):
    """
      This class define an interactive plot area with which the user
      can dynamically set the data by clicking within the plotting area
    """
    
    def __init__(self, parent=None,Uxdata=[0],Uydata=[0],Thxdata=[0],Thydata=[0]):
        self.fig = Figure()
        super(InteractivePlot, self).__init__(self.fig)
        self.setParent(parent)
        
        # Create axes for the figure
        # ax1 is for Uranium
        # ax2 is for Thorium
        self.ax1 = self.figure.add_subplot(111)
        self.ax1.minorticks_on()
        self.ax1.tick_params(axis='y', labelcolor='red')
        self.ax1.grid(visible=True,which='major', color='0.8', linestyle='-')
        self.ax1.grid(visible=True,which='minor', color='0.5', linestyle='--')
        self.ax2 = self.ax1.twinx()
        self.ax2.tick_params(axis='y', labelcolor='blue')
        self.ax2.grid(visible=True,which='major', color='0.8', linestyle='-')
        self.ax2.grid(visible=True,which='minor', color='0.5', linestyle='--')
        
        # Store points
        self.listPoints = []
        # self.mpl_connect('button_press_event', self.plot_draggable_point)
        
        # Create the plot
        self.bottomLeftX = 0
        self.bottomLeftY = 0
        self.topRightX = 0
        self.topRightY = 0
        self.x = np.array(
            [
                self.bottomLeftX,
                self.bottomLeftX,
                self.topRightX,
                self.topRightX,
                self.bottomLeftX,
                ])
        self.y = np.array(
            [
                self.bottomLeftY,
                self.topRightY,
                self.topRightY,
                self.bottomLeftY,
                self.bottomLeftY,])
        (self.myPlot,) = self.ax1.plot(self.x, self.y)

        self.ax1.set_xlim((0,1))
        self.ax1.set_ylim((0,100))
        self.ax2.set_ylim((0,200))
        
        # Plot the concentrations
        self.Uxdata = np.array(Uxdata)
        self.Uydata = np.array(Uydata)
        self.Thxdata = np.array(Thxdata)
        self.Thydata = np.array(Thydata)

        self.ax1.set_xlabel('edge <- Grain distance -> core')
        self.ax1.set_ylabel('U (ppm)', color = 'red')
        self.ax2.set_ylabel('Th (ppm)', color = 'blue')
        
        # Create draggable points
        self.create_draggable_points()
        self.fig.canvas.mpl_connect('axes_enter_event', self.on_enter)
        self.fig.canvas.mpl_connect('key_press_event', self.on_press)
        
        self.plotSnap = 0.1
    
    #------------------------------------------------------
    def on_enter(self, _):
        # When the mouse pointer enter in the plot area
        self.fig.canvas.setFocus()
        
    #------------------------------------------------------
    def create_draggable_points(self):
          self.listPoints.append(DraggablePoint(self, True, self.Uxdata, self.Uydata, 5, color='r', ax = self.ax1))
          self.listPoints.append(DraggablePoint(self, False, self.Thxdata, self.Thydata, 5, color='b', ax = self.ax2))
          # self.newPointLabel = "point2"
          # self.listPoints.append(DraggablePoint(self, True, self.xdata[1], self.ydata[1], 5))
          
    #-----------------------------------------------------
    def on_press(self, event):
        # When press one of the following keyuu
        if event.key == 'u':
            self.ax1.set_visible(True)
            self.ax2.set_visible(False)
            self.listPoints[0].enable = True
            self.listPoints[1].enable = False
            self.draw()
            
        if event.key == 't':
            self.ax1.set_visible(False)
            self.ax2.set_visible(True)
            self.ax2.get_xaxis().set_visible(True)
            self.ax2.set_xlabel('Grain distance')
            self.listPoints[0].enable = False
            self.listPoints[1].enable = True
            self.draw()
            
        if event.key == 'r':
            self.ax1.set_visible(True)
            self.ax2.set_visible(True)
            self.listPoints[0].enable = False
            self.listPoints[1].enable = False
            self.draw()
            
    #------------------------------------------------------
    def updateFigure(self):
       self.draw()
       
       
#############################################################
class MessageBox(QWidget):
    """ 
    This class set a question message box and center it in the interface
    
    """
    def __init__(self, title='PecubeGUI message box', text1='message'):
        super().__init__()
        self.title = title
        self.left = 10
        self.top = 10
        self.height = 400
        self.width = 600
        self.text1 = text1
        self.initUI()
    
    #----------------------------------------------------------
    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        IconPecube = os.path.join(conf.IconPath,"Pecube_icon.ico")
        self.setWindowTitle('New project')
        self.setWindowIcon(QIcon(IconPecube))
        #▬Force the window to open at the center of the screen
        qtRectangle = self.frameGeometry()
        centerPoint = QDesktopWidget().availableGeometry().center()
        qtRectangle.moveCenter(centerPoint)
        self.move(qtRectangle.topLeft())
        
        self.QMessage = QMessageBox.question(self,'PecubeGUI message box', self.text1,
                                              QMessageBox.Yes | QMessageBox.Cancel | QMessageBox.No, QMessageBox.No)
        
        
#############################################################
class WindowProject(QWidget, object):
    """ 
    This class ask the user a name for the new Pecube project.
    """
    def __init__(self, parent, ProjName):
        super().__init__()
        IconNewInput = os.path.join(conf.IconPath,"New_Input.png")
        self.setWindowTitle('New project')
        self.parent = parent
        self.setWindowIcon(QIcon(IconNewInput))
        self.setGeometry(800, 400, 200, 100)
        #force the window to open at the center of the screen     
        qtRectangle = self.frameGeometry()
        centerPoint = QDesktopWidget().availableGeometry().center()
        qtRectangle.moveCenter(centerPoint)
        self.move(qtRectangle.topLeft())
        self.pname = ProjName.ProjectName
        self.ProjectName = self.getText()
    
    #----------------------------------------------------------
    def getText(self):
        """ The dialog window that asks the user a name for the Pecube project
        is handled here
        """
        text, okPressed = QInputDialog.getText(
            self, "New project", "Enter a project Name (5 characters):", QLineEdit.Normal, self.pname)
        if len(text) != 5 and okPressed == True:
            QErrorMessage(self.parent).showMessage('Incorrect project name size, please keep 5 characters.')
        elif okPressed and text != '': 
            self.ProjectName = str(text)
            return self.ProjectName
        else:
            #The user quit the window or cancel
            self.ProjectName = 'thisisempty'
            return self.ProjectName