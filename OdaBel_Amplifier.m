function varargout = OdaBel_Amplifier(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @OdaBel_Amplifier_OpeningFcn, ...
                   'gui_OutputFcn',  @OdaBel_Amplifier_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before OdaBel_Amplifier is made visible.
function OdaBel_Amplifier_OpeningFcn(hObject, eventdata, handles,varargin)
    handles.output = hObject;
    guidata(hObject, handles);
    
    % ---------- STATUS UPDATE ------------ %
    set(handles.status,'String', 'Status: ');
    
    axes(handles.logo_plot);
        imshow(imread('logo.jpg'), []);

    axes(handles.plot_input);
        title('Input wave');         
        xlabel('Time [s]'); 
        ylabel('Amplitude [V]');
        grid on;    

    axes(handles.plot_tri_input);
        title('Comparator inputs');   
        xlabel('Time [s]'); 
        ylabel('Amplitude [V]');
        grid on;

    axes(handles.plot_pwm);
        title('Comparator output'); 
        xlabel('Time [s]'); 
        ylabel('Amplitude [V]');
        grid on;

    axes(handles.plot_output_time);
        title('Amplifier output');  
        xlabel('Time [s]'); 
        ylabel('Amplitude [V]');
        grid on;

    axes(handles.plot_pwm_filter);
        title('Comparator output FT & Filter response FT');   
        xlabel('Frequency [Hz]'); 
        ylabel('Amplitude [dB]');
        grid on;

    axes(handles.plot_output_freq);
        title('Amplifier output FT');        
        xlabel('Frequency [Hz]'); 
        ylabel('Amplitude [dB]');
        grid on;

    numSlider = (get(handles.slider,'Value'));
    set(handles.textSlider,'String',round(numSlider,1));
    

% --- Executes on slider movement.
function slider_Callback(hObject, eventdata, handles)
    numSlider = (get(handles.slider,'Value'));
    set(handles.textSlider,'String',round(numSlider,1));
    
% --- Executes on button press in run_button.
function run_button_Callback(hObject, eventdata, handles)

    %% -------------------------------------------------------------------    
    % AXES RESET
    % --------------------------------------------------------------------

    axes(handles.plot_input);       cla reset     
    axes(handles.plot_tri_input);   cla reset   
    axes(handles.plot_pwm);         cla reset   
    axes(handles.plot_output_time); cla reset   
    axes(handles.plot_pwm_filter);  cla reset   
    axes(handles.plot_output_freq); cla reset

    % ---------- STATUS UPDATE -------------- %
    set(handles.status,'String', 'Status: 0%');


    %% -------------------------------------------------------------------    
    % VARIABLES
    % --------------------------------------------------------------------

    fsin    = str2double(get(handles.freq_input, 'String'));
    Asin    = str2double(get(handles.a_input,    'String'));
    ftri    = str2double(get(handles.freq_triang,'String'));
    G       = str2double(get(handles.gain,       'String'));
    R       = str2double(get(handles.resistencia,'String'));
    L       = str2double(get(handles.bobina,     'String')) * 1e-6 ;
    C       = str2double(get(handles.capacitat,  'String')) * 1e-9;
    Nperiod = 250;  % Factor to obtain Fs
    Fs      = ftri * Nperiod;
    T       = 0.25; % Time duration
    Ts      = 1/Fs;
    margin  = Asin * 0.2; % Plot amplitude margin for better visualization
    % Logo colors
    purple  = [0.328 0.209 0.296];
    red     = [0.852 0.289 0.296]; 


    % Quality factor
    Q       = num2str( round( R*sqrt(C/L), 2) );
    set(handles.Q,'String', Q);
    set(handles.Q,'BackgroundColor',[0.5 0.5 0.5]);

    % Cutoff frequency
    calcul_fc = num2str( round((sqrt(1/(L*C*4*pi^2)))/1000 ,2) );
    set(handles.fc,'String',calcul_fc);
    set(handles.fc,'BackgroundColor',[0.5 0.5 0.5]);


    %% -------------------------------------------------------------------    
    % INPUT SIGNAL
    % --------------------------------------------------------------------

    tsin    = 0:Ts:T;
    sinWave = Asin*sin(2*pi*fsin*tsin);

    axes(handles.plot_input);
        plot(tsin,sinWave,'LineWidth',1,'Color',red);
        signal_xlim = get(handles.slider,'Value');
        xlim([0 signal_xlim*(1/fsin)]);
        ylim([-Asin-margin Asin+margin]);
        title('Input wave');        
        grid on; 
        xlabel('Time [s]'); 
        ylabel('Amplitude [V]');

    % ---------- STATUS UPDATE --------------- %
    set(handles.status,'String', 'Status: 10%');


    %% -------------------------------------------------------------------    
    % TRIANGULAR WAVE
    % --------------------------------------------------------------------

    % Create a larger rectangular wave in order to avoid the non-periodic 
    % part of the next convolution
    rect    = 0:Ts:2*T;
    j       = 1;
    for i = 1:length(rect)
        if(j <= Nperiod/2)
            rect(i) = 1;
        else
            rect(i) = -1;
        end
        if (j == Nperiod)
            j = 1;
        else 
            j = j + 1;
        end
    end

    % ---------- STATUS UPDATE --------------- %
    set(handles.status,'String', 'Status: 20%');

    unRect  = rect(1: Nperiod);
    triWave = conv(rect,unRect);
    % Cut the initial and final part of the convolution
    % Ensuring that length(triWave) = length(tsin)
    triWave = triWave(Nperiod*2 : Nperiod*2+length(tsin)-1);
    triWave = Asin*triWave/max(triWave);


    %% -------------------------------------------------------------------    
    % COMPARATOR
    % --------------------------------------------------------------------

    axes(handles.plot_tri_input);
        plot(tsin,triWave,'LineWidth',1,'Color',purple); hold on
        plot(tsin,sinWave,'LineWidth',1,'Color',red);
        xlim([0 signal_xlim*(1/fsin)]);
        ylim([-Asin-margin Asin+margin]);
        title('Comparator inputs');  
        xlabel('Time [s]'); 
        ylabel('Amplitude [V]');
        grid on; 

    % ---------- STATUS UPDATE --------------- %
    set(handles.status,'String', 'Status: 30%');
    
    outComp = 0:Ts:T;
    for i = 1:length(sinWave)
       if (triWave(i) < sinWave(i))
           outComp(i) = -Asin;
       else
           outComp(i) = Asin;
       end    
    end
    % PWM signal

    % Gain
    outGain = outComp * G;

    axes(handles.plot_pwm);
        plot(tsin,outComp,'LineWidth',1,'Color',red);
        xlim([0 signal_xlim*(1/fsin)]);
        ylim([(min(outComp)-margin) (max(outComp)+margin)]);
        title('Comparator output');
        xlabel('Time [s]'); 
        ylabel('Amplitude [V]');
        grid on;

    % ---------- STATUS UPDATE --------------- %
    set(handles.status,'String', 'Status: 40%');


    % PWM FT
    N       = 2^24;
    FTpwm   = (2/N)*abs(fft(outGain,N));
    FTpwm   = FTpwm + 1e-50; % Avoid zeros in order to apply the logarithm
    FTpwmdB = 20*log10(FTpwm);

    % ---------- STATUS UPDATE --------------- %
    set(handles.status,'String', 'Status: 50%');


    %% -------------------------------------------------------------------    
    % LOW-PASS FILTER
    % --------------------------------------------------------------------

    % H(z) numerator
    b = [Ts^2*R 2*Ts^2*R Ts^2*R];
    % H(z) denominator
    x = (R*C*L*4)+(2*Ts*L)+(Ts^2*R);    % Independent term
    y = (2*Ts^2*R)-(8*R*C*L);           % z^-1 term
    z = (R*C*L*4)-(2*Ts*L)+(Ts^2*R);    % z^-2 term
    a = [x y z];


    % Impulsional response
    delta       = zeros(1,length(tsin));
    delta(1)    = 1;
    respImp     = filter(b,a,delta);
    FTrespImp   = abs(fft(respImp,N));
    FTrespImpdB = 20*log10(FTrespImp);

    % ---------- STATUS UPDATE --------------- %
    set(handles.status,'String', 'Status: 60%');

    f = Fs*(1/N:1/N:1);
    axes(handles.plot_pwm_filter);     
        semilogx(f,FTpwmdB,    'LineWidth',1,'Color',red); hold on;  
        semilogx(f,FTrespImpdB,'LineWidth',1,'Color',purple);   
        xlim([20 ftri*4+ftri/2]);           
        ylim([-80 max(FTpwmdB+10)]);
        title('Comparator output FT & Filter response FT');   
        xlabel('Frequency [Hz]'); 
        ylabel('Amplitude [dB]');
        grid on;


    %% -------------------------------------------------------------------    
    % OUTPUT SIGNAL
    % --------------------------------------------------------------------

    Vo = filter(b,a,outGain);

    % ---------- STATUS UPDATE --------------- %
    set(handles.status,'String', 'Status: 70%');

    axes(handles.plot_output_time);
        plot(tsin,Vo,'LineWidth',1,'Color',red);
        xlim([0 signal_xlim*(1/fsin)]);
        ylim([(-max(Vo)-margin*G)  (max(Vo)+margin*G)]);
        title('Amplifier output');
        xlabel('Time [s]'); 
        ylabel('Amplitude [V]');
        grid on;

    % ---------- STATUS UPDATE --------------- %
    set(handles.status,'String', 'Status: 80%');

    FTVo    = (2/N)*abs(fft(Vo,N));
    FTVodB  = 20*log10(FTVo);

    % ---------- STATUS UPDATE --------------- %
    set(handles.status,'String', 'Status: 90%');

    axes(handles.plot_output_freq);   
        semilogx(f,FTVodB,'LineWidth',1,'Color',red);
        xlim([20 ftri*4+ftri/2]);           
        ylim([-80 max(FTVodB+10)]);
        title('Amplifier output FT');       
        xlabel('Frequency [Hz]'); 
        ylabel('Amplitude [dB]');
        grid on; 

    % ---------- STATUS UPDATE ---------------- %
    set(handles.status,'String', 'Status: 100%');











%% -------------------------------------------------------------------    
% Functions created by default for each element of the interface:
% --------------------------------------------------------------------

% Executes during object creation, after setting all properties.
function freq_triang_CreateFcn(hObject, eventdata, handles)
% hObject    handle to freq_triang (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function a_triang_Callback(hObject, eventdata, handles)
% hObject    handle to a_triang (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a_triang as text
%str2double(get(hObject,'String')) returns contents of a_triang as a double

% --- Executes during object creation, after setting all properties.
function a_triang_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a_triang (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function freq_input_Callback(hObject, eventdata, handles)
    
function a_input_Callback(hObject, eventdata, handles)

function freq_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to freq_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function gain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to freq_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function a_input_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a_input (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function freq_triang_Callback(hObject, eventdata, handles)
    
function varargout = OdaBel_Amplifier_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;
    
function gain_Callback(hObject, eventdata, handles)

function resistencia_Callback(hObject, eventdata, handles)
% hObject    handle to resistencia (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of resistencia as text
%str2double(get(hObject,'String')) returns contents of resistencia as a double

% --- Executes during object creation, after setting all properties.
function resistencia_CreateFcn(hObject, eventdata, handles)
% hObject    handle to resistencia (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function bobina_Callback(hObject, eventdata, handles)
% hObject    handle to bobina (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bobina as text
%str2double(get(hObject,'String')) returns contents of bobina as a double

% --- Executes during object creation, after setting all properties.
function bobina_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bobina (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function capacitat_Callback(hObject, eventdata, handles)
% hObject    handle to capacitat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of capacitat as text
% str2double(get(hObject,'String')) returns contents of capacitat as a double

% --- Executes during object creation, after setting all properties.
function capacitat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to capacitat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function Q_Callback(hObject, eventdata, handles)
% hObject    handle to Q (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Q as text
%        str2double(get(hObject,'String')) returns contents of Q as a double

% --- Executes during object creation, after setting all properties.
function Q_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Q (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function fc_Callback(hObject, eventdata, handles)
% hObject    handle to fc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fc as text
%        str2double(get(hObject,'String')) returns contents of fc as a double

% --- Executes during object creation, after setting all properties.
function fc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
