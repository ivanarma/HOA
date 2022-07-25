"""
@author: Ivan ARMAND
"""
import itertools #to concatenate lists of lists
import numpy as np #dynamic lists
import matplotlib.pyplot as plt #to plot curves
from scipy.optimize import least_squares #to solve non-linear equations
from scipy.optimize import curve_fit
import tkinter as tk #to create a window with frames and canvas
from tkinter import filedialog #to search a file
from tkinter import ttk #to use some widgets
import pandas as pd #to read xlsx and clipboard
import xlsxwriter
from matplotlib.figure import Figure #to draw figures on tkinter frame
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import warnings
import webbrowser
warnings.filterwarnings("ignore")

####### Non-linear fitting functions ##################### 

P_2=["Hill-Deboer","Fowler-Guggenheim","Langmuir","Freundlich","Temkin","dubinin-Radushkevich","Flory-Huggins","Hill","Halsey","Harkin-Jura","Jovanovic","Elovich","Kiselev"]
P_3=["Redlich-Peterson","Sips","Toth","Koble-Carrigan","Kahn","Radke-Prausniiz","Langmuir-Freundlich","Jossens"]
P_4=["Fritz-Schlunder","Baudu","Marczewski-Jaroniec"]
Equations=[['affine', 'y=ax+b'], ["Henry's", 'y=ax'], ['Langmuir', 'y=x/(a+x/b)'], ['Freundlich', 'y=ln(y)=ln(a)+ln(x)/b'], ['Temkin', 'y=a+b*ln(x)'], ['Halsey', 'y=(b-ln(x))/a'], ['Harkin-Jura', 'y=abs(1/(b/a-(1/a)log10(x)))**0.5'], ['Jovanovic', 'ln(y)=ln(a)-b*x'], ['Redlich-Peterson', 'y=a*x/(1+b*x^c)'], ['Sips', 'y=(a*x^c)/(1-b*x^c)'], ['Toth', 'y=a*x/((1+b*x^c)^(1/c))'], ['Koble-Carrigan', 'y=1/((1/(a*x^b))+c/a)'], ['Kahn', 'y=(a*b*x)/((1+b*x)*c)'], ['Radke-Prausniiz', 'y=a*b*x/((1+b*x)^c)'], ['Langmuir-Freundlich', 'y=(a*(b*x)^c)/(1+(b*x)^c)'], ['Fritz-Schlunder', 'y=a*x^b/(1+c*x^d)'], ['Marczewski-Jaroniec', 'y=a*(((b*x)^c)/(1+(b*x)^c))^d']]
analysis_method="virial"
virial_number_of_a_parameters=5
virial_number_of_b_parameters=2

def virial_function(X,*constants):#ln(y)=ln(x)+sum(ai*x^i)/T+sum(bi*x^i)
    if isinstance (constants[0], list ):
        constants=constants[0]
    lena,lenb=5,2
    x,T=X
    s=np.log(x)
    for i in range(lena):
        s+=(constants[i]*x**i)/T
    for i in range(lena,lena+lenb):
        s+=(constants[i]*x**(i-lena))
    return s

def curve_fit_virial(x,T,Ydatas):
    """we assume that T has the same length as x and Ydatas are the results of f(x,T)"""
    p0=[1,1,1,1,1,1,1]
    return curve_fit(virial_function,(x,T),Ydatas,p0)[0]

def G_modelfunction(x0,X,Y,fit="Langmuir-Freundlich"):
    """contain the model equation for different fitting methods, initial conditions are needed and initial experimental arrays too. If we return 0 the fit will be bad so we won't keep it"""
    a=x0[0]
    b=x0[1]
    c=x0[2]
    d=x0[3]
    if fit=="affine":
        return a+b*X-Y
    ###1 parameter fit
    if fit=="Henry's": #y=ax
        return a*X-Y
    
    ###2 parameters fits
    if fit =="Langmuir": #y=X/(a+X/b)
        return X/(a+X/b)-Y
    if fit =="Freundlich": #ln(y)=ln(a)+ln(x)/b
        return np.exp(np.log(a)+np.log(X)/b)-Y
    if fit =="Temkin": #y=a+b*ln(x)
        return a+b*np.log(X)-Y
    if fit=="Halsey":#y=(b-ln(x))/a
        return (b-np.log(X))/a-Y
    if fit=="Harkin-Jura":#y=abs(1/(b/a-(1/a)log10(x)))**0.5
        return abs(1/(b/a-(np.log10(X)/a)))**0.5-Y
    if fit=="Jovanovic": #ln(y)=ln(a)-b*x
        return np.exp(np.log(a)-b*X)-Y
    
    ###3 parameters fits
    if fit=="Redlich-Peterson": #y=a*x/(1+b*x^c)
        return a*X/(1+b*X**c)-Y
    if fit =="Sips": #y=(a*x**c)/(1-b*x^c)
        return (a*X**c)/(1+b*X**c)-Y
    if fit =="Toth": #y=a*x/((1+b*x^c)^(1/c))
        return a*b*X/((1+(b*X)**c)**(1/c))-Y
    if fit=="Koble-Carrigan": #y=1/((1/(a*x^b))+c/a)
        return 1/((1/(a*X**b))+c/a)-Y
    if fit=="Kahn": #y=(a*b*x)/((1+b*x)*c)
        return (a*b*X)/((1-b*X)*c)-Y
    if fit=="Radke-Prausniiz":#y=a*b*x/((1+b*x)^c)
        return (a*b*X)/((1+b*X)**c)-Y
    if fit=="Langmuir-Freundlich": #y=(a*(b*x)^c)/(1+(b*x)^c)
        return (a*(b*X)**c)/(1+(b*X)**c)-Y
    
    ###4 parameters fits
    if fit=="Fritz-Schlunder": #y=a*x^b/(1+c*x^d)
        return (a*X**b)/(1+c*X**d)-Y
    if fit=="Marczewski-Jaroniec": #y=a*(((b*x)^c)/(1+(b*x)^c))^d
        return a*(((b*X)**c)/(1+(b*X)**c))**d-Y

def G_Q(Xf,Gx,fit="Langmuir-Freundlich"): #Y=G_Q(pp,lsqnonlin(GP,GQ).x) Gx-->x0 & Xf-->X
    """this method return the results of the fitting methods, Gx are the constants that were found, we want to have the Y's for the Xf's values"""
    a=Gx[0]
    b=Gx[1]
    c=Gx[2]
    d=Gx[3]
    
    if fit=="affine":
        return a+b*Xf
    
    ###1 parameter fit    
    if fit=="Henry's":
        return a*Xf
    
    ###2 parameters fits
    if fit =="Langmuir":
        return Xf/(a+Xf/b)
    if fit=="Freundlich":
        return 10**(np.log(a)+np.log(Xf)/b)
    if fit =="Temkin":
        return a+b*np.log(Xf)
    if fit=="Halsey":
        return (b-np.log(Xf))/a
    if fit=="Harkin-Jura":
        return abs(1/(b/a-(np.log10(Xf)/a)))**0.5
    if fit=="Jovanovic": #ln(y)=ln(a)-b*x
        return np.exp(np.log(a)-b*Xf)
    
    ###3 parameters fits
    if fit=="Redlich-Peterson":
        return a*Xf/(1+b*Xf**c)
    if fit=="Sips":
        return (a*Xf**c)/(1+b*Xf**c)
    if fit =="Toth":
        return a*b*Xf/((1+(b*Xf)**c)**(1/c))
    if fit=="Koble-Carrigan":
        return 1/((1/(a*Xf**b))+c/a)
    if fit=="Kahn":
        return (a*b*Xf)/((1-b*Xf)*c)
    if fit=="Radke-Prausniiz":
        return (a*b*Xf)/((1+b*Xf)**c)
    if fit=="Langmuir-Freundlich":
        return (a*(b*Xf)**c)/(1+(b*Xf)**c)
    
    ###4 parameters fits
    if fit=="Fritz-Schlunder":
        return (a*Xf**b)/(1+c*Xf**d)
    if fit=="Marczewski-Jaroniec":
        return a*(((b*Xf)**c)/(1+(b*Xf)**c))**d

def lsqnonlin(X,Y,fit="Langmuir-Freundlich",x0=np.array([1,1,1,1])):
    """return the constants of the model function using the initial conditions and the experimental arrays"""
    return least_squares(G_modelfunction, x0, args=(X, Y,fit)) 

def fit(Xf,X,Y,fit="Langmuir-Freundlich",x0=np.array([1,1,1,1])):
    """fit the array X,Y. For an Xf array given Yf is returned, f(Xf)=Yf"""
    Yf=G_Q(Xf,lsqnonlin(X,Y,fit,x0).x,fit)
    return Yf

def find_best_fits(Columns):
    """Compare the mean of RMSE for the 3curves, For each fitting method"""
    RMSEy=[]
    x0=np.array([x0a.get(),x0b.get(),x0c.get(),x0d.get()])
    for method in Fits[1:]:
        if method in P_2:
            pass
        elif method in P_3:
            pass
        elif method in P_4:
            pass
        else:
            pass
        RMSEx=[]
        for j in range(2,len(Columns),4):
            X,Y=Columns[j],Columns[j+1]
            Xf=X
            Yf=fit(Xf, X, Y,method,x0)
            RMSEx.append(RMSE(Yf-Y))
        RMSEy.append(np.mean(RMSEx))
    index_min=np.argmin(RMSEy)
    return Fits[1:][index_min]
                
########################################################

####### ERROR ANALYSIS FUNCTIONS ##################### 
    
def RMSE(A):
    """Root mean square errors"""
    return np.sqrt(np.mean(A**2))

def chi_square(Observed,Calculated):
    s=0
    for i in range(len(Observed)):
        s+=((Calculated[i]-Observed[i])**2)/Observed[i]
    return s

def Coefficient_of_nondetermination(Observed,Calculated):
    s1,s2=0,0
    qmobs=np.mean(Observed)
    for i in range(len(Observed)):
        s1+=(Calculated[i]-Observed[i])**2
        s2+=(Calculated[i]-Observed[i])**2+(Calculated[i]-qmobs)**2
    return 1-s1/s2

def ERRSQ(Observed,Calculated):
    """Sum square of errors"""
    s=0
    for i in range(len(Observed)):
        s+=(Observed[i]-Calculated[i])**2
    return s

def EASB(Observed,Calculated):
    """Sum of absolute errors"""
    s=0
    for i in range(len(Observed)):
        s+=abs(Observed[i]-Calculated[i])
    return s

def HYBRID(Observed,Calculated,p):
    """Hybrid fractional error function, p=number of parameters of the fit function"""
    s=0
    for i in range(len(Observed)):
        s+=((Calculated[i]-Observed[i])**2)/Observed[i]
    return 100*s/(len(Observed)-p)

def MPSD(Observed,Calculated,p):
    """Marquardt’s Percent Standard Deviation"""
    s=0
    for i in range(len(Observed)):
        s+=((Calculated[i]-Observed[i])**2)/Observed[i]
    return (s/(len(Observed)-p))**0.5

def ARE(Observed,Calculated):
    """Average Relative Error"""
    s=0
    for i in range(len(Observed)):
        s+=(Calculated[i]-Observed[i])/Observed[i]
    return 100*s/len(Observed)

########################################################
    
####### file reading functions #######################

def isfloat(num):
    """detect if num is an integer or a float"""
    try:
        float(num)
        return True
    except ValueError:
        return False

def split_column(C):
    """split column C into two columns, we suppose that there are at least one empty row between the two columns"""
    C1,C2=[],[]
    C=np.append(C," ")
    i=0
    first_column=True
    while(i<len(C)):
        if isfloat(str(C[i]).replace(',','.')) and str(C[i])!="nan": #if it is a number
            if first_column:
                C1.append(float(str(C[i]).replace(',','.')))
            else:
                C2.append(float(str(C[i]).replace(',','.')))
        else:
            if i>5: #meaning we are sure that it was not the header
                first_column=False
        i+=1
    return [C1,C2]

def open_folder(): #we want to know where to put the results
    """get access to a folder path"""
    return filedialog.askdirectory()

def open_file():
    """get access to a file path"""
    global paste_button
    file_path = filedialog.askopenfilename(filetypes = (("excel files", "*.xlsx"),("All files", "*.*")))
    if file_path[-5:]==".xlsx":
        A=pd.read_excel(file_path, index_col=None, header=None)
        Col=[split_column(A.values[0:,i]) for i in range(int(Number_of_isotherms.get())*2)] #hence we have the N*2 columns we want as numpy arrays
        load_datas(Col)

def paste_clipboard():
    """Two kind or inputs are accepted, the first is one column containing adsorption and desorption below with at least one empty row between the two columns. The second is just values of adsorption"""
    global datas
    global paste_button
    clipboard=pd.read_clipboard(sep='\t')
    datasx=clipboard.values[:,0]
    datasy=clipboard.values[:,1]
    Colx=split_column(datasx)
    Coly=split_column(datasy)
    datas.append(Colx)
    datas.append(Coly)
    onetime=True
    for i in range(len(Temperatures)-1):
        if paste_button.cget('text')=="paste datas from the "+str(Temperatures[i])+"°K isotherm" and onetime:
            paste_button.config(text="paste datas from the "+str(Temperatures[i+1])+"°K isotherm")
            onetime=False
    if paste_button.cget('text')=="paste datas from the "+str(Temperatures[-1])+"°K isotherm" and onetime:
        paste_button.destroy()
        load_datas(datas)

def load_datas(Col):
    global analysis_method
    """according to the 6columns given-read-, we sort them into a Columns variable. We suppose that Colx[0] means adsorption, and colx[1] means desorption, colx[1] are unused in this program, we just see it because it simplify inputs"""
    datas_frame.destroy()
    analysis_method=method_analysis_selector.get()
    method_selector_frame.destroy()
    if analysis_method=="clausius":
        fitting_frame.pack()
    if analysis_method=="virial":
        plot_resultsx.destroy()
    plot_results.pack()
    plot_frame.pack()
    Columns=[]
    for i in range(0,len(Col),2):
        Columns.append(Col[i][1])
        Columns.append(Col[i+1][1])
        Columns.append(Col[i][0])
        Columns.append(Col[i+1][0])
    Columns=[np.array(i) for i in Columns]
    plot_button = tk.Button(master = plot_frame, command = lambda: plot_curves(Columns),height = 2, width = 10,text = "Plot")
    unplot_button = tk.Button(master = plot_frame, command =clear_plot,height = 2, width = 10,text = "Unplot")
    save_button = tk.Button(master = plot_frame, command =lambda: save_results(Columns),height = 2, width = 10,text = "Save results")
    write_errors_button = tk.Button(master = plot_frame, command =lambda: draw_errors_analysis(Columns),height = 2, width = 10,text = "errors_array")
    draw_array(Columns)
    plot_button.pack(side="left")
    unplot_button.pack(side="left")
    save_button.pack(side="left")
    write_errors_button.pack(side="left")
    paste_button.destroy()

def write_results(Cols,folder_path,best_fit):
    df=pd.DataFrame()
    writer=pd.ExcelWriter(folder_path+"/"+best_fit+" HOA datas.xlsx",engine="xlsxwriter")
    df.to_excel(writer,sheet_name="Sheet1")
    workbook=writer.book
    worksheet=writer.sheets["Sheet1"]
    Nf,HOA=Cols
    Nf,HOA=np.array(Nf),np.array(HOA)
    worksheet.write(0,0,"Amount adsorbed n [mmol/g]")
    worksheet.write(0,1,"-ΔHads [kJ/mol]")
    worksheet.write(0,2,"made by ARMAND Ivan")
    worksheet.write(0,3,"https://github.com/ivanarma/HOA")
    for i in range(len(Nf)):
        worksheet.write(i+1,0,"="+str(Nf[i]))
        worksheet.write(i+1,1,"="+str(HOA[i]))
    chart = workbook.add_chart({'type': 'line'})
    chart.add_series({
    'name': 'HOA',
    'categories': '=Sheet1!$A$2:$A$'+str(len(Nf)+1),
    'values':     '=Sheet1!$B$2:$B$'+str(len(Nf)+1),
})
    chart.set_x_axis({'name': 'Amount adsorbed n [mmol/g]'})
    chart.set_y_axis({'name': '-ΔHads [kJ/mol]','major_gridlines': {'visible': False}})
    worksheet.insert_chart('D2', chart)
    workbook.close()

def save_results(Columns):
    """create a png file for every graph in the selected folder. Be aware that if a file with the same name exists, it replace the file, so there is a need of one folder for one dataset"""
    if analysis_method=="clausius":
        best_fit=fit_selector.get()
        if fit_selector.get()=="automatic":
            best_fit=find_best_fits(Columns)
        folder_path=open_folder()
        Nf,HOA=save_graph5(Columns, folder_path,best_fit)
        save_graph4(Columns, folder_path,best_fit)
        save_graph3(Columns, folder_path,best_fit)
        save_graph2(Columns, folder_path,best_fit)
        save_graph1(Columns, folder_path,best_fit)
        write_results((Nf,HOA),folder_path,best_fit)
    if analysis_method=="virial":
        folder_path=open_folder()
        Nf,HOA=save_graph_virial(Columns, folder_path)
        write_results((Nf,HOA),folder_path,"virial")

def findX_from_Y(y, xdata, ydata):
    """find the x in xdata corresponding to the y in ydata, we suppose that there is only one solution, but it stays an approximation as we just have arrays, not functions. We suppose that the curve is descending. We use dichotomy"""
    xl,xr=0,len(ydata) #l in xl stands for left and r in xr for right
    x=xr//2
    for i in range(int(np.log2(len(ydata)))):
        yval=ydata[x]-y
        if yval<0:
            xl=x
            x+=(xr-xl)//2
        else:
            xr=x
            x-=(xr-xl)//2    
    return xdata[x]

####### Graphical-related functions #####################

def number_to_string(value,nsd):
    """convert the number to string with nsd as number of significant digits"""
    return str(("{:."+str(nsd)+"f}").format(round(value,nsd)))

def create_table(method,Values,p):
    """give the errors table for a fitting method"""
    global tableframe
    popup=tk.Tk()
    popup.title(method+" error table")
    Functions=["RMSE","R²","X²","ERRSQ","EASB","HYBRID","MPSD","ARE","a","b","c","d"]
    if p!=4:
        Functions=Functions[:-4+p]
    tableframe=tk.Frame(popup)
    tableframe.pack(fill = tk.BOTH,expand=tk.YES)
    scrollv = tk.Scrollbar(tableframe,orient="vertical")
    scrollv.pack(side="right", fill="y")
    scrollh = tk.Scrollbar(tableframe,orient='horizontal')
    scrollh.pack(side="bottom", fill="x")
    table=ttk.Treeview(tableframe,yscrollcommand=scrollv.set, xscrollcommand = scrollh.set)
    table.column("#0", width=0,  stretch="no")
    scrollv.config(command=table.yview)
    scrollh.config(command=table.xview)
    for tab in Equations:
            if tab[0]==method:
                tk.Label(tableframe, text = "method : "+method+", "+tab[1]).pack(side="top")
    if method=="virial":
        Functions=["RMSE","R²","X²","ERRSQ","EASB","HYBRID","MPSD","ARE"]+["a"+str(i) for i in range(1,virial_number_of_a_parameters+1)]+["b"+str(i) for i in range(1,virial_number_of_b_parameters+1)]
        tk.Label(tableframe, text = "virial analysis").pack(side="top")
        fig = Figure(figsize=(7, 1), dpi=100)
        wx = fig.add_subplot(111)
        wx.text(0.2, 0.6, r'$ ln(p) = ln(n) + \frac{1}{T} \sum_{i=0}^{4} a_{i+1} n^{i} + \sum_{i=0}^{1} b_{i+1} n^{i}$', fontsize = 20)
        canvas = FigureCanvasTkAgg(fig, master=tableframe)
        canvas.get_tk_widget().pack(side="top", fill=tk.BOTH, expand=1)
        canvas._tkcanvas.pack(side="top", fill=tk.BOTH, expand=1)
        # Set the visibility of the Canvas figure
        wx.get_xaxis().set_visible(False)
        wx.get_yaxis().set_visible(False)
        
        table["columns"]=("col1","col2")
        table.column("col1",anchor="center",width=80)
        table.heading("col1",text="ln(y)=ln(x)+sum(ai*x^i)/T+sum(bi*x^i)",anchor="center")
        table.column("col2",anchor="center", width=80)
        table.heading("col2",text="value",anchor="center")
        for f in range(0,len(Functions)):
            val=[Functions[f],Values[f]]

            table.insert(parent='',index='end',iid=f,text='',values=tuple(val))
    else:
        table["columns"]=tuple([method]+[str(T) for T in Temperatures])
        table.column(method,anchor="center", width=80)
        table.heading(method,text=method,anchor="center")
        for T in range(len(Temperatures)):
            table.column(str(Temperatures[T]),anchor="center", width=150)
            table.heading(str(Temperatures[T]),text=str(Temperatures[T]),anchor="center")
        for f in range(0,len(Functions)):
            val=[Functions[f]]
            for T in range(len(Temperatures)):
                val.append(Values[T][f])
            table.insert(parent='',index='end',iid=f,text='',values=tuple(val))
    table.pack(fill = tk.BOTH,expand=tk.YES)
    menubar = tk.Menu(popup)
    popup.config(menu=menubar)
    file_menu = tk.Menu(menubar,tearoff=False)
    file_menu.add_command(label='Exit',command=popup.destroy)
    menubar.add_cascade(label="File",menu=file_menu,underline=0)
    popup.mainloop()

def draw_errors_analysis(Columns):
    """when you click error-array button """
    method=fit_selector.get()
    graph.delete("err")
    Values=[]
    if fit_selector.get()=="automatic":
        method=find_best_fits(Columns)
    p=1
    if method in P_2:
        p=2
    if method in P_3:
        p=3
    if method in P_4:
        p=4
    if analysis_method=="virial":
        p=virial_number_of_a_parameters+virial_number_of_b_parameters
        method="virial"
        XX,YY=[],[]
        for j in range(2,len(Columns),4):
            GP,GQ=Columns[j],Columns[j+1]
            GP,GQ=GQ,np.log(GP*1000)
            XX.append(GP)
            YY.append(GQ)
        lenx=np.array([len(xx) for xx in XX])
        TT=[]
        for T in range(len(Temperatures)):
            for i in range(lenx[T]):
                TT.append(int(Temperatures[T]))
        XX=np.array(list(itertools.chain.from_iterable(XX)))
        YY=np.array(list(itertools.chain.from_iterable(YY)))
        constants=curve_fit_virial(XX, TT, YY)
        print(constants)
        Xf=XX
        a1,a2,a3,a4,a5,b1,b2=constants[0],constants[1],constants[2],constants[3],constants[4],constants[5],constants[6]
        Yf=virial_function((XX,TT), a1,a2,a3,a4,a5,b1,b2)
        Y=YY
        rmse_t=number_to_string(RMSE(Yf-Y), nsd)
        R2_t=number_to_string(Coefficient_of_nondetermination(Y, Yf), nsd)
        X2_t=number_to_string(chi_square(Y, Yf), nsd)
        errsq_t=number_to_string(ERRSQ(Y, Yf), nsd)
        easb_t=number_to_string(EASB(Y,Yf), nsd)
        hybrid_t=number_to_string(HYBRID(Y,Yf,p), nsd)
        mpsd_t=number_to_string(MPSD(Y,Yf,p), nsd)
        are_t=number_to_string(ARE(Y,Yf), nsd)
        Values=[rmse_t,R2_t,X2_t,errsq_t,easb_t,hybrid_t,mpsd_t,are_t,constants[0],constants[1],constants[2],constants[3],constants[4],constants[5],constants[6]]
    else:
        for j in range(2,len(Columns),4):
            X,Y=Columns[j],Columns[j+1]
            Xf=X
            Yf=fit(Xf, X, Y,method)
            constants=lsqnonlin(X, Y,method).x
            Val=[]
            for k in range(len(constants[:p])):
                c=number_to_string(constants[k], nsd)
                Val.append(c)
            rmse_t=number_to_string(RMSE(Yf-Y), nsd)
            R2_t=number_to_string(Coefficient_of_nondetermination(Y, Yf), nsd)
            X2_t=number_to_string(chi_square(Y, Yf), nsd)
            errsq_t=number_to_string(ERRSQ(Y, Yf), nsd)
            easb_t=number_to_string(EASB(Y,Yf), nsd)
            hybrid_t=number_to_string(HYBRID(Y,Yf,p), nsd)
            mpsd_t=number_to_string(MPSD(Y,Yf,p), nsd)
            are_t=number_to_string(ARE(Y,Yf), nsd)
            Values.append([rmse_t,R2_t,X2_t,errsq_t,easb_t,hybrid_t,mpsd_t,are_t]+Val)
    create_table(method,Values,p)

def draw_array(Columns):
    """draw the table containing all the experimental datas entered"""
    px=26
    col_name=[]
    for i in range(0,len(Columns),4):
        col_name.append("Xdes"+str(i//4+1))
        col_name.append("Ydes"+str(i//4+1))
        col_name.append("Xads"+str(i//4+1))
        col_name.append("Yads"+str(i//4+1))
    for i in range(len(Columns)):
        for j in range(len(Columns[i])):
            graph.create_text(int(szx*(100+3*px/4+i*px*2)),int(szy*(120+10+j*px)),text="{:.4f}".format(round(Columns[i][j], 4)),font=('Arial',fontsize,'normal'))
            graph.create_rectangle(int(szx*(96+i*px*2)),int(szy*(120+j*px)),int(szx*(96+2*px*(i+1))),int(szy*(120+px*(j+1))),)
        graph.create_text(int(szx*(100+3*px/4+i*px*2)),int(szy*(120+10-px)),text=col_name[i],font=('Arial',fontsize,'normal'))

def plot_graph1(Columns,fig=Figure(figsize = (5, 5),dpi = 100)):
    """Adsorption+Desorption curves"""
    plotx = fig.add_subplot(111)
    plotx.set_title("Isotherms")
    for j in range(0,len(Columns)-1,2): #first graph 2min40s
        plotx.plot(Columns[j],Columns[j+1],linestyle=' ',marker='o')
    legende=[]
    for T in Temperatures:
        legende.append("Des"+str(T)+"°K")
        legende.append("Ads"+str(T)+"°K")
    plotx.legend(legende)
    plotx.set_xlabel("Pressure (kPa)")
    plotx.set_ylabel("Amount adsorbed (mmol/g)")
    return fig

def save_graph1(Columns,foldername,method="Langmuir-Freundlich"):
    """Adsorption+Desorption curves"""
    plt.clf()
    plt.title("Isotherms")
    for j in range(0,len(Columns)-1,2): #first graph 2min40s
        plt.plot(Columns[j],Columns[j+1],linestyle=' ',marker='o')
    legende=[]
    for T in Temperatures:
        legende.append("Des"+str(T)+"°K")
        legende.append("Ads"+str(T)+"°K")
    plt.legend(legende)
    plt.xlabel("Pressure (kPa)")
    plt.ylabel("Amount adsorbed (mmol/g)")
    plt.savefig(foldername+"/adsorption+desorption curves.png")

def plot_graph2(Columns,fig=Figure(figsize = (5, 5),dpi = 100)):
    """adsorption curves"""
    plotx = fig.add_subplot(111)
    plotx.set_title("Only adsorption isotherms")
    for j in range(2,len(Columns),4): #second graph 3min2s
        plotx.plot(Columns[j],Columns[j+1],linestyle=' ',marker='o') #GP,GQ
    legende=[]
    for T in Temperatures:
        legende.append("Ads"+str(T)+"°K")
    plotx.legend(legende)
    plotx.set_xlabel("Pressure (kPa)")
    plotx.set_ylabel("Amount adsorbed (mmol/g)")
    return fig

def save_graph2(Columns,foldername,method="Langmuir-Freundlich"):
    """adsorption curves"""
    plt.clf()
    plt.title("Only adsorption isotherms")
    for j in range(2,len(Columns),4): #second graph 3min2s
        plt.plot(Columns[j],Columns[j+1],linestyle=' ',marker='o') #GP,GQ
    legende=[]
    for T in Temperatures:
        legende.append("Ads"+str(T)+"°K")
    plt.legend(legende)
    plt.xlabel("Pressure (kPa)")
    plt.ylabel("Amount adsorbed (mmol/g)")
    plt.savefig(foldername+"/adsorption curves.png")
    
def plot_graph3(Columns,method="Langmuir-Freundlich",fig=Figure(figsize = (5,5),dpi = 100,)):
    """Non-linear fit curves"""
    fig.clear()
    plotx = fig.add_subplot(111)
    plotx.set_title("Fit of adsorption curves, fit : "+method)
    N=np.arange(float(interva.get()),float(intervb.get()),(float(intervb.get())-float(interva.get()))/int(Number_of_amount_adsorbed_points.get()))
    legende=[]
    pp=np.linspace(0,1,int(Number_of_fitted_curve_points.get()))#number of points of the fitted curve
    XX=[]
    for j in range(2,len(Columns),4):
        GP,GQ=Columns[j],Columns[j+1]
        Y=fit(pp, GP, GQ,method)
        X=[]
        for i in N:
            X.append(findX_from_Y(i, pp, Y))
        XX.append(X)
        plotx.plot(pp,Y)
        legende.append("fit curve "+str(Temperatures[j//4])+"°K")
    for j in range(2,len(Columns),4): #second graph curves another time
        plotx.plot(Columns[j],Columns[j+1],linestyle=' ',marker='o') #GP,GQ
    for t in Temperatures:
        legende.append("Ads"+str(t)+"°K")
    plotx.legend(legende)  
    plotx.set_xlabel("Pressure (kPa)")
    plotx.set_ylabel("Amount adsorbed (mmol/g)")
    return fig

def save_graph3(Columns,foldername,method="Langmuir-Freundlich"):
    plt.clf()
    plt.title("Fit of adsorption curves, fit : "+method)
    N=np.arange(float(interva.get()),float(intervb.get()),(float(intervb.get())-float(interva.get()))/int(Number_of_amount_adsorbed_points.get()))#4min20
    legende=[]
    pp=np.linspace(0,0.1,int(Number_of_fitted_curve_points.get()))#number of points of the fitted curve
    XX=[]
    for j in range(2,len(Columns),4):
        GP,GQ=Columns[j],Columns[j+1]
        Y=fit(pp, GP, GQ,method)
        X=[]
        for i in N:    
            X.append(findX_from_Y(i, pp, Y))
        XX.append(X)
        plt.plot(pp,Y)
        legende.append("fit curve "+str(Temperatures[j//4])+"°K")
    for j in range(2,len(Columns),4): #second graph curves another time
        plt.plot(Columns[j],Columns[j+1],linestyle=' ',marker='o') #GP,GQ
    for t in Temperatures:
        legende.append("Ads"+str(t)+"°K")
    plt.legend(legende)  
    plt.xlabel("Pressure (kPa)")
    plt.ylabel("Amount adsorbed (mmol/g)")
    plt.savefig(foldername+"/Isotherm fit "+method+" method.png")

def plot_graph4(Columns,method="Langmuir-Freundlich",fig=Figure(figsize = (5, 5),dpi = 100)):
    """loading curves"""
    fig.clear()
    plotx = fig.add_subplot(111)
    plotx.set_xlabel("1/T (°K)")
    plotx.set_ylabel("ln p")
    plotx.set_title("Curves for each n (amount adsorbed or loading) in mmol/g")
    XX=[]
    pp=np.linspace(0,1,int(Number_of_fitted_curve_points.get()))#number of points of the fitted curve
    N=np.arange(float(interva.get()),float(intervb.get()),(float(intervb.get())-float(interva.get()))/int(Number_of_amount_adsorbed_points.get()))#4min20
    for j in range(2,len(Columns),4):
        GP,GQ=Columns[j],Columns[j+1]
        Y=fit(pp, GP, GQ,method)
        X=[]
        for i in N:
            X.append(findX_from_Y(i,pp,Y))
        XX.append(X)
    legende=[]
    for n in range(len(N)):
        x=[]
        y=[]
        for t in range(len(Temperatures)):
            if XX[t][n]!=None:
                x.append(XX[t][n])
                y.append(1/Temperatures[t])
        legende.append("n = "+str(int(N[n]*int(Number_of_amount_adsorbed_points.get()))/int(Number_of_amount_adsorbed_points.get())))
        x,y=np.array(y),np.array(x)
        y=np.log(y)
        plotx.plot(x,y)
    plotx.legend(legende)
    return fig

def save_graph4(Columns,foldername,method="Langmuir-Freundlich"):
    """loading curves"""
    plt.clf()
    plt.xlabel("1/T (°K)")
    plt.ylabel("ln p")
    plt.title("Curves for each n (amount adsorbed or loading) in mmol/g")
    XX=[]
    pp=np.linspace(0,1,int(Number_of_fitted_curve_points.get()))#number of points of the fitted curve
    N=np.arange(float(interva.get()),float(intervb.get()),(float(intervb.get())-float(interva.get()))/int(Number_of_amount_adsorbed_points.get()))#4min20
    for j in range(2,len(Columns),4):
        GP,GQ=Columns[j],Columns[j+1]
        Y=fit(pp, GP, GQ,method)
        X=[]
        for i in N:
            X.append(findX_from_Y(i,pp,Y))
        XX.append(X)
    legende=[]
    for n in range(len(N)):
        x=[]
        y=[]
        for t in range(len(Temperatures)):
            if XX[t][n]!=None:
                x.append(XX[t][n])
                y.append(1/Temperatures[t])
        legende.append("n = "+str(int(N[n]*int(Number_of_amount_adsorbed_points.get()))/int(Number_of_amount_adsorbed_points.get())))
        x,y=np.array(y),np.array(x)
        y=np.log(y)
        plt.plot(x,y)
    plt.legend(legende)
    plt.savefig(foldername+"/ln P vs T-1 "+method+" method.png")

def plot_graph5(Columns,method="Langmuir-Freundlich",fig=Figure(figsize = (5, 5),dpi = 100)):
    """HOA Curve"""
    fig.clear()
    plotx = fig.add_subplot(111)
    pp=np.linspace(0,1,int(Number_of_fitted_curve_points.get()))#number of points of the fitted curve
    N=np.arange(float(interva.get()),float(intervb.get()),(float(intervb.get())-float(interva.get()))/int(Number_of_amount_adsorbed_points.get()))#4min20
    plotx.set_title("Isosteric heat/enthalpy of adsorption (HOA), fit : "+method)
    plotx.set_xlabel("Amount adsorbed n [mmol/g]")
    plotx.set_ylabel("-ΔHads [kJ/mol]")
    XX=[]
    for j in range(2,len(Columns),4):
        GP,GQ=Columns[j],Columns[j+1]
        Y=fit(pp, GP, GQ,method)
        X=[]
        for i in N:
            X.append(findX_from_Y(i,pp,Y))
        XX.append(X)
    HOA=[]
    Nf=[]
    for n in range(len(N)):
        x=[]
        y=[]
        h=0
        for t in range(len(Temperatures)):
            if XX[t][n]!=None:
                h+=1
                x.append(XX[t][n])
                y.append(1/Temperatures[t])
        if h>1:
            Nf.append(N[n])
            x,y=np.array(y),np.array(x)
            y=np.log(y)
            pp=np.array([1/Temperatures[-1],1/Temperatures[0]])
            Y=fit(pp,x,y,"affine")
            m=(Y[1]-Y[0])/(pp[1]-pp[0])
            HOA.append(-m*8.314/1000)      
    plotx.plot(Nf,HOA)
    plotx.legend(["HOA"])
    return fig

def save_graph5(Columns,foldername,method="Langmuir-Freundlich"):
    """HOA Curve"""
    plt.clf()
    pp=np.linspace(0,1,int(Number_of_fitted_curve_points.get()))#number of points of the fitted curve
    N=np.arange(float(interva.get()),float(intervb.get()),(float(intervb.get())-float(interva.get()))/int(Number_of_amount_adsorbed_points.get()))#4min20
    plt.title("Isosteric heat/enthalpy of adsorption (HOA), fit : "+method)
    plt.xlabel("Amount adsorbed n [mmol/g]")
    plt.ylabel("-ΔHads [kJ/mol]")
    XX=[]
    for j in range(2,len(Columns),4):
        GP,GQ=Columns[j],Columns[j+1]
        Y=fit(pp, GP, GQ,method)
        X=[]
        for i in N:
            X.append(findX_from_Y(i,pp,Y))
        XX.append(X)
    HOA=[]
    Nf=[]
    for n in range(len(N)):
        x=[]
        y=[]
        h=0
        for t in range(len(Temperatures)):
            if XX[t][n]!=None:
                h+=1
                x.append(XX[t][n])
                y.append(1/Temperatures[t])
        if h>1:
            Nf.append(N[n])
            x,y=np.array(y),np.array(x)
            y=np.log(y)
            pp=np.array([1/Temperatures[-1],1/Temperatures[0]])
            Y=fit(pp,x,y,"affine")
            m=(Y[1]-Y[0])/(pp[1]-pp[0])
            HOA.append(-m*8.314/1000)
    plt.plot(Nf,HOA)
    plt.legend(["HOA"])
    plt.savefig(foldername+"/HOA "+method+" method.png")
    return Nf,HOA

def plot_graph_virial(Columns,fig=Figure(figsize=(5,5),dpi=100)):
    global virial_constants
    fig.clear()
    plotx = fig.add_subplot(111)
    plotx.set_title("Fit of adsorption curves, virial analysis "+r'$ ln(p) = ln(n) + \frac{1}{T} \sum_{i=0}^{4} a_{i+1} n^{i} + \sum_{i=0}^{1} b_{i+1} n^{i}$')
    N=np.arange(float(interva.get()),float(intervb.get()),(float(intervb.get())-float(interva.get()))/int(Number_of_amount_adsorbed_points.get()))
    XX=[]
    YY=[]
    for j in range(2,len(Columns),4):
        GP,GQ=Columns[j],Columns[j+1]
        GP,GQ=GQ,np.log(GP*1000)
        XX.append(GP)
        YY.append(GQ)
    lenx=np.array([len(xx) for xx in XX])
    TT=[]
    for T in range(len(Temperatures)):
        for i in range(lenx[T]):
            TT.append(int(Temperatures[T]))
    XX=np.array(list(itertools.chain.from_iterable(XX)))
    YY=np.array(list(itertools.chain.from_iterable(YY)))
    constants=curve_fit_virial(XX, TT, YY)
    HOA=0
    for i in range(5):
        HOA+=8.314*constants[i]*N**i
    HOA=HOA/1000
    plotx.plot(N,-HOA)
    plotx.legend(["HOA virial"])  
    plotx.set_xlabel("Amount adsorbed (mmol/g)")
    plotx.set_ylabel("-ΔHads [kJ/mol]")
    return fig

def save_graph_virial(Columns,foldername):
    plt.clf()
    plt.title("Heat of adsorption curve, virial analysis")
    Nf=np.arange(float(interva.get()),float(intervb.get()),(float(intervb.get())-float(interva.get()))/int(Number_of_amount_adsorbed_points.get()))
    XX=[]
    YY=[]
    for j in range(2,len(Columns),4):
        GP,GQ=Columns[j],Columns[j+1]
        GP,GQ=GQ,np.log(GP*1000)
        XX.append(GP)
        YY.append(GQ)
    lenx=np.array([len(xx) for xx in XX])
    TT=[]
    for T in range(len(Temperatures)):
        for i in range(lenx[T]):
            TT.append(int(Temperatures[T]))
    XX=np.array(list(itertools.chain.from_iterable(XX)))
    YY=np.array(list(itertools.chain.from_iterable(YY)))
    constants=curve_fit_virial(XX, TT, YY)
    HOA=0
    for i in range(5):
        HOA+=8.314*constants[i]*Nf**i
    HOA=-HOA/1000
    plt.plot(Nf,HOA)
    plt.legend(["HOA virial"])  
    plt.xlabel("Amount adsorbed (mmol/g)")
    plt.ylabel("-ΔHads [kJ/mol]")
    plt.savefig(foldername+"/HOA_virial.png")
    return Nf,HOA
    
def draw_plot(fig,x=0,y=0):
    """will add the defined fig in the pan from right to left"""
    canvas=FigureCanvasTkAgg(fig,master=pan)
    canvas.get_tk_widget().pack(side=tk.RIGHT,fill=tk.BOTH,expand=tk.YES)
    canvas.draw()

def clear_plot():
    """unplot"""
    for child in pan.winfo_children():
        if child!=graph:
            child.destroy()
    
def plot_curves(Columns):
    """here we find the fitting method, then we plot what we want"""
    clear_plot()
    if analysis_method=="clausius":
        best_fit=fit_selector.get()
        if fit_selector.get()=="automatic":
            best_fit=find_best_fits(Columns)
        draw_plot(plot_graph3(Columns,best_fit))
        draw_plot(plot_graph5(Columns,best_fit))
    if analysis_method=="virial":
        draw_plot(plot_graph_virial(Columns))

def mousewheel(event):
    """move the text if you use the mousewheel, if you can't see everything due to widgets or graphs"""
    graph.yview_scroll(-int(event.delta/120*szy), "units") 

def setTemperature(T):
    global Temperatures
    Temperatures=T
    paste_button['text']="paste datas from the "+str(T[0])+"°K isotherm"
    datas_frame.destroy()

def asking_temperatures_of_curves():
    set_temperatures_button['text']="set temperatures (°K)"
    entryNum4.destroy()
    Ti=[]
    for i in range(int(Number_of_isotherms.get())):
        tk.Label(datas_frame, text = "T"+str(i)).pack(side="left")
        Ti.append(tk.StringVar())
        entryNum4i=tk.Entry(datas_frame, textvariable=Ti[i])
        if i<len(Temperatures):
            entryNum4i.insert(0,Temperatures[i])
        entryNum4i.pack(side="left")
    set_temperatures_button.configure(command=lambda : setTemperature([int(t.get()) for t in Ti]))

def setting_analysis_method(name):
    global analysis_method
    analysis_method=name

def how_it_works():
    popup=tk.Tk()
    popup.title('HOA calculus about page')
    popup.geometry('600x400+50+50')   
    menubar = tk.Menu(popup)
    popup.config(menu=menubar)
    file_menu = tk.Menu(menubar,tearoff=False)
    file_menu.add_command(label='Exit',command=popup.destroy)
    menubar.add_cascade(label="File",menu=file_menu,underline=0)
    pan=tk.Frame(popup)
    pan.pack(expand=tk.YES,fill=tk.BOTH)
    Line1="Select analysis method : clausius-clapeyron or virial analysis, default : clausius-clapeyron"
    Line2="Select number of isotherms you possess for calculus, then push the 'set_curves_number&temperatures' button, default : 3isotherms"
    Line3="Then enter the temperatures your isotherms are about in kelvin. Please the numbers should be sorted in ascending order, then push the 'set_temperatures' button, default : 273,283,293 °K"
    Line4="There are 2 ways to enter datas, "
    Line5="You have to copy in your clipboard two excel column at a time containing adsorption isotherm (left column = relative pressure, right column = quantity adsorbed (mmol/g)) , if so, please push 'paste button'"
    Line6="You have an excel file containing the datas obtained from the software MicroActive when you push 'copy data' of an isotherm plot and you assemble everything from left to right without space"
    Line7="If everything is correct, according to the method used the interface will change a little "
    Line8="For clausius clapeyron : you can choose what fitting equations the program will use, default : automatic"
    Line9="You can choose the initial values for the fitting methods (it is a bit useless)"
    Line10="You can choose the scale of the HOA plot (don't choose less than 0 or equal to 0 for clausius-clapeyron) - default [0.3,2.3]"
    Line11="You can choose the number of points of the fitting curves, and of the horizontal axis of HOA. The more the longer the program will take but the more precise the results will be"
    Line12="You can click to the'errors_array' button to have the results of the quality of the selected fit for each isotherm fit in a table"
    Line13="You can save the results. Create a folder, go to this folder and then select it. 5 png files will be generated and 1 excel file "
    Line14="For virial analysis : 'errors_array' we display the 7 parameters, for HOA calculus only the 5 a parameters are important"
    Line15="7 parameters are used ln(y)=ln(x)+sum(ai*x^i)/T+sum(bi*x^i) with a1,a2,a3,a4,a5 and b1,b2, two variable x and T"
    Line16="Save results will generate 1png file and one excel corresponding of the HOA curve obtained"
    Line17="23 June 2022, program made by Ivan ARMAND during internship at Ångström laboratory."
    url="https://github.com/ivanarma?tab=repositories"
    Nums=["1.","2.1.","2.2.","3.","3.A","3.B","4.","4.A.1","4.A.2","4.A.3","4.A.4","4.A.5","4.A.6","4.B.1","4.B.2","4.B.3","date of the build/author"]
    Text=[Line1,Line2,Line3,Line4,Line5,Line6,Line7,Line8,Line9,Line10,Line11,Line12,Line13,Line14,Line15,Line16,Line17]
    lent=[len(l) for l in Text]
    m=max(lent)+500
    for line in range(len(Text)):
        while len(Text[line])<=m:
            Text[line]+=" "
        lframe=tk.Frame(pan)
        tk.Label(lframe, text = Nums[line],fg='red',anchor="w").pack(side="left")
        tk.Label(lframe, text = Text[line],anchor="w").pack(side="right")
        lframe.pack()
    lframe=tk.Frame(pan)
    tk.Label(lframe,text="Please check my github page for more informations : ",foreground="DarkOrchid3",font= ('Aerial 18')).pack(side="left")
    link=tk.Label(lframe,text=url,cursor= "hand2", foreground= "dodger blue",font= ('Aerial 18'))
    link.bind("<Button-1>", lambda e:open_url(url))
    link.pack(side="right")
    lframe.pack()
    popup.mainloop()

def open_url(url):
   webbrowser.open_new_tab(url)
########################################################    

###creating a window###
root = tk.Tk()
root.title('HOA calculus')
root.geometry('600x400+50+50')
#add panel
pan = tk.Frame(root) #
pan.pack(expand = tk.YES,fill = tk.BOTH) #panel is filling the window
graph = tk.Canvas(pan,bg='gray') #add a canvas
graph.pack(fill = tk.BOTH,expand=tk.YES,side="bottom") #canvas is filling the panel
graph.bind_all("<MouseWheel>", mousewheel)
#######################


###adding a menu###
menubar = tk.Menu(root)
root.config(menu=menubar)
file_menu = tk.Menu(menubar,tearoff=False)
file_menu.add_command(label='open_folder',command=open_folder)
file_menu.add_command(label='open xlsx file',command=open_file)
file_menu.add_separator()
file_menu.add_command(label='Exit',command=root.destroy)
about_menu = tk.Menu(menubar,tearoff=False)
about_menu.add_command(label='how it works ?',command=how_it_works)
menubar.add_cascade(label="File",menu=file_menu,underline=0)
menubar.add_cascade(label="about the program",menu=about_menu,underline=0)
####################

###adding widgets###
nsd=7 #number of significant digits for every number displayed
Temperatures=[273,283,293]#default values


###adding frames###
plot_frame=tk.Frame(graph)
fitting_frame=tk.Frame(graph)
datas_frame=tk.Frame(graph)
plot_results=tk.Frame(graph)
method_selector_frame=tk.Frame(pan)
datas_frame.pack()
###

###method selector

tk.Label(method_selector_frame,text="analysis method").pack(side="left")
method_analysis_selector=ttk.Combobox(method_selector_frame,values=["clausius","virial"])
method_analysis_selector.current(0)
method_analysis_selector.pack(side="right")
method_selector_frame.pack()
###
###fit selector-->fitting frame
tk.Label(fitting_frame,text="Fitting equation").pack()
Fits=["automatic"]+[eq[0] for eq in Equations[1:]] #we don't consider affine as an isotherm method
fit_selector = ttk.Combobox(fitting_frame,values=Fits)
fit_selector.current(0)
fit_selector.pack()
###

###fitting related widgets-->fitting frame
x0a = tk.StringVar()
x0b = tk.StringVar()
x0c = tk.StringVar()
x0d = tk.StringVar()
tk.Label(fitting_frame, text = "a0=").pack(side="left")
entryNum4a=tk.Entry(fitting_frame, textvariable=x0a)
entryNum4a.insert(0,1)
entryNum4a.pack(side="left")
tk.Label(fitting_frame, text = "; b0=").pack(side="left")
entryNum4b=tk.Entry(fitting_frame, textvariable=x0b)
entryNum4b.insert(0,1)
entryNum4b.pack(side="left")
tk.Label(fitting_frame, text = "; c0=").pack(side="left")
entryNum4c=tk.Entry(fitting_frame, textvariable=x0c)
entryNum4c.insert(0,1)
entryNum4c.pack(side="left")
tk.Label(fitting_frame, text = "; d0=").pack(side="left")
entryNum4d=tk.Entry(fitting_frame, textvariable=x0d)
entryNum4d.insert(0,1)
entryNum4d.pack()
###
###plots related widgets-->plot_results frame

plot_resultsx=tk.Frame(plot_results)
plot_resultsx.pack(side="bottom")
Number_of_fitted_curve_points = tk.StringVar()
tk.Label(plot_resultsx, text = "Number of points of the fitted curves").pack(side="left")
entryNum1 = tk.Entry(plot_resultsx, textvariable=Number_of_fitted_curve_points)
entryNum1.insert(0,10000)
entryNum1.pack(side="left")

pan3=tk.Frame(plot_results)
pan3.pack(side="bottom")
Number_of_amount_adsorbed_points = tk.StringVar()
tk.Label(pan3, text = "Number of points of the 'amount adsorbed' axis").pack(side="left")
entryNum2 = tk.Entry(pan3, textvariable=Number_of_amount_adsorbed_points)
entryNum2.insert(0,100)
entryNum2.pack(side="left")

pan4=tk.Frame(plot_results,borderwidth=3, relief="ridge")
pan4.pack(side="bottom")
interva = tk.StringVar()
intervb = tk.StringVar()
tk.Label(pan4, text = "Interval of 'amount adsorbed' axis [").pack(side="left")
entryNum3a = tk.Entry(pan4, textvariable=interva)
entryNum3a.insert(0,0.3)
entryNum3a.pack(side="left")
tk.Label(pan4, text = ";").pack(side="left")
entryNum3b = tk.Entry(pan4, textvariable=intervb)
entryNum3b.insert(0,2.3)
entryNum3b.pack(side="left")
tk.Label(pan4, text = "]").pack(side="left")
paste_button = tk.Button(master = graph, command = paste_clipboard,height = 2, width = 30,text = "paste datas from the "+str(Temperatures[0])+"°K isotherm",borderwidth=3, relief="solid")
paste_button.pack(side="top")
###

###datas related widgets -->datas_frame
Number_of_isotherms = tk.StringVar()
entryNum4=tk.Entry(datas_frame, textvariable=Number_of_isotherms)
entryNum4.insert(0,3)
entryNum4.pack(side="right")
set_temperatures_button=tk.Button(master = datas_frame, command = asking_temperatures_of_curves,text = "set_curves_number&temperatures",borderwidth=3, relief="solid")
set_temperatures_button.pack(side="left")

###
datas=[]
szx,szy=root.winfo_screenwidth()/1536,root.winfo_screenheight()/864
fontsize=int(10/szx)-1
graph.yview_scroll(int(100*szy), "units") 
tableframe=None
####################
root.mainloop()