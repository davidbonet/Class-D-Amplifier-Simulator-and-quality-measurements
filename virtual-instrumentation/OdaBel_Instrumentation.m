function varargout = OdaBel_Instrumentation(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @OdaBel_Instrumentation_OpeningFcn, ...
                   'gui_OutputFcn',  @OdaBel_Instrumentation_OutputFcn, ...
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

% --- Executes just before OdaBel_Instrumentation is made visible.
function OdaBel_Instrumentation_OpeningFcn(hObject, eventdata, handles, varargin)
    handles.output = hObject;
    guidata(hObject, handles);

    % ---------- STATUS UPDATE ------------ %
    set(handles.status,'String', 'Status: ');
    
    axes(handles.logoplot);
        imshow(imread('logo.jpg'), []);
    
    axes(handles.THD_time);
        title('Amplifier output (THD)');   
        xlabel('Time [s]'); 
        ylabel('Amplitude');
        grid on; 
    
    axes(handles.THD_freq);
        title('Amplifier output FT (THD)'); 
        xlabel('Frequency [Hz]'); 
        ylabel('Amplitude');
        grid on;   
    
    axes(handles.freq_response);
        title('Amplifier frequency response');
        xlabel('Frequency [Hz]'); 
        ylabel('Amplitude [dB]');
        grid on;

    axes(handles.out_ampli_fr);
        title('Amplifier output (FR)');   
        xlabel('Time [s]'); 
        ylabel('Amplitude');
        grid on; 
        
    set(handles.freq_margin,'String',get(handles.sliderfreq,'Value'));
    

% --- Executes on slider movement.
function sliderfreq_Callback(hObject, eventdata, handles)
    set(handles.freq_margin,'String',get(handles.sliderfreq,'Value'));    
    

function varargout = OdaBel_Instrumentation_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;



% --- Executes on button press in thdbutton.
function thdbutton_Callback(hObject, eventdata, handles)

    %% -------------------------------------------------------------------    
    % RESET
    % --------------------------------------------------------------------    

    % ---------- STATUS UPDATE -------------- %
    set(handles.status,'String', 'Status: 0%');

    axes(handles.THD_time);     cla reset     
    axes(handles.THD_freq);     cla reset
    set (handles.THDresult,'String','');


    %% -------------------------------------------------------------------    
    % VARIABLES
    % --------------------------------------------------------------------    

    Fs      = str2double(get(handles.editFs,'String'));
    Asin    = str2double(get(handles.editA_sin,'String'));
    f_in    = str2double(get(handles.editf_sin,'String'));
    N       = 2^22;
    Ts      = 1/Fs;
    T       = 2;    % Time duration
    t       = 0:Ts:T;
    Nbits   = 16;
    NChan   = 1;
    margin  = 0.2; % Plot amplitude margin for better visualization
    % Logo colors
    red     = [0.852 0.289 0.296];
    purple  = [0.328 0.209 0.296];


    %% -------------------------------------------------------------------    
    % INPUT SIGNAL
    % --------------------------------------------------------------------    

    signalIn    = Asin*sin(2*pi*f_in*t);

    % Send and receive amplifier's output signal (DA/AD)
    player      = audioplayer(signalIn,Fs,Nbits);
    data        = audiorecorder(Fs,Nbits,NChan);
    play(player);
    recordblocking(data,T);
    signalOut   = getaudiodata(data);

    if(max(abs(signalOut))>0.95)
         warndlg('Turn down the microphone volume','Saturated signal');
         return
    end

    
    %% -------------------------------------------------------------------    
    % PROCESS RECEIVED SIGNAL
    % --------------------------------------------------------------------    

    signalOut   = transpose(signalOut);
    inici       = 0.5*Fs;   % Discard first samples
    final       = (T-0.5)*Fs;
    signalOut   = signalOut(inici : final);      

    t = 0:Ts:1;
    axes(handles.THD_time);
        plot(t,signalOut,'Color',red,'LineWidth',1); 
        title('Amplifier output (THD)');   
        xlabel('Time [s]'); 
        ylabel('Amplitude');
        xlim([0 8*(1/f_in)]);
        ylim([-margin-max(signalOut) max(signalOut)+margin]); 
        grid on;

    % ---------- STATUS UPDATE --------------- %
    set(handles.status,'String', 'Status: 25%');

    % Use Hamming window
    load('ventanaHamming');
    L = 2^13;
    window = transpose(hammingWindowToolbox);
    outWindow  = signalOut(L:L*2-1).*window;

    % ---------- STATUS UPDATE --------------- %
    set(handles.status,'String', 'Status: 50%');

    outF    = abs(fft(outWindow,N));
    F_fund  = f_in/Fs;
    k_fund  = round(F_fund*N);  % First harmonic (F = k/N)
    AmpFund = outF(k_fund);     % Amplitude fundamental frequency
    outF    = 1/AmpFund*outF;   % Normalize

    % Frequency plot
    f = Fs*(1/N:1/N:1);
    axes(handles.THD_freq);
        semilogx(f,outF,'Color',red,'LineWidth',1);
        title('Amplifier output FT (THD)'); 
        xlabel('Frequency [Hz]'); 
        ylabel('Amplitude');
        xlim([20 Fs/2]);
        grid on
        hold on

    % ---------- STATUS UPDATE --------------- %
    set(handles.status,'String', 'Status: 75%');


    %% -------------------------------------------------------------------    
    % CALCULATE THD
    % --------------------------------------------------------------------    
    i = 2;
    harmonics = 0;
    k = k_fund*2;
    k_last = round(0.5*N); % DFT sample of Fs/2

    % Go through all the harmonics until we reach Fs/2
    while k < k_last
        harmonics = harmonics + (outF(k))^2;

        scatter(f(k),outF(k),10,'MarkerEdgeColor',purple);
        hold on

        i = i+1;
        k = k_fund*i;
    end

    harmonics = sqrt(harmonics); % THD Formula (Normalized amplitude)
    THD = harmonics*100;
    set(handles.THDresult,'Value',THD);
    percent = '%';
    str = sprintf('THD: %.3f%s',THD,percent);
    set(handles.THDresult,'String',str);
    set(handles.THDresult, 'BackgroundColor', [0.5 0.5 0.5]);
    set(handles.THDresult, 'ForegroundColor', [1 1 1]);

    % ---------- STATUS UPDATE ---------------- %
    set(handles.status,'String', 'Status: 100%');







% --- Executes on button press in freqbutton.
function freqbutton_Callback(hObject, eventdata, handles)
    %% -------------------------------------------------------------------    
    % RESET
    % --------------------------------------------------------------------    

    % ---------- STATUS UPDATE -------------- %
    set(handles.status,'String', 'Status: 0%');

    axes(handles.freq_response);    cla reset     
    axes(handles.out_ampli_fr);     cla reset    

    amptype = get(handles.type,'Value'); % 1 = Half-Bridge, 2 = Full-Bridge
    check_sat = get(handles.checksaturation,'Value');
   
    THD = get(handles.THDresult,'Value');
    if  THD == 0.0 && check_sat == 1
        warndlg('Calculate THD to ensure that the signal is not being sent saturated','Check THD');
        return
    elseif THD > 1.0 && check_sat == 1
        warndlg('THD higher than expected, check the volumes of the PC','Check volumes');
        return
    end


    %% -------------------------------------------------------------------    
    % VARIABLES
    % --------------------------------------------------------------------

    Fs          = str2double(get(handles.editFs,'String'));
    Nfreq       = str2double(get(handles.numfreq, 'String'));
    duration    = 1;         % Records duration time in seconds;
    Num_Chan    = 2;         % 1 Means Mono. 2 Means Stereo;
    Nbits       = 16;        % Number of bits of A/D D/A converters;
    t           = (1/Fs):(1/Fs):duration;
    vectGaindB  = zeros(1,Nfreq);
    gain        = zeros(1,Nfreq);
    freq        = logspace(log10(10),log10((Fs/2)*0.9),Nfreq);
    marginfreq  = get(handles.sliderfreq,'Value');
    margintime  = 0.01;
    red         = [0.852 0.289 0.296];
    data_ampli  = audiorecorder(Fs, Nbits,Num_Chan);


    %% -------------------------------------------------------------------    
    % CALCULATE FREQUENCY RESPONSE
    % --------------------------------------------------------------------

    for i = 1:Nfreq
        axes(handles.out_ampli_fr);
        input = 0.8*sin(2*pi*freq(i)*t); % Generate diferent sinusoids
        
        % Send and receive amplifier's output signal (DA/AD)
        in_ampli = audioplayer(input,Fs,16);      
        play(in_ampli);       
        recordblocking(data_ampli, duration);

        output  = getaudiodata(data_ampli);
        out_L   = output(:,1);
        out_R   = output(:,2);

            plot(t,out_R,'LineWidth',1,'Color',red);
            xlim([0.9 1]);
            ylim([-margintime+min(out_R) margintime+max(out_R)]);
            title('Amplifier output (FR)'); 
            xlabel('Time [s]'); 
            ylabel('Amplitude');
            grid on;

        if amptype == 1 % Half-Bridge
            power_L = out_L'*out_L;
            power_R = out_R'*out_R;
            gain(i) = power_R/power_L;
            vectGaindB(i) = 10*log10(gain(i));
        else % Full-Bridge
            output = out_L - out_R; % Differential
            power = output'*output;
            gain(i) = power;
            vectGaindB(i) = 10*log10(gain(i));
        end

        axes(handles.freq_response);
            semilogx(freq(i),vectGaindB(i),'-o','Color',red,'MarkerSize',4,...
            'MarkerFaceColor',red); hold on;
            ylim ([min(vectGaindB)-marginfreq max(vectGaindB)+marginfreq]); 
            xlim([10 Fs/2*0.9]);
            ylabel('Gain (dB)');
            xlabel('Frequency (Hz)');
            title('Frequency response');
            grid on;

        status  = round((i/Nfreq)*100);
        percent = '%';
        str2 = sprintf('Status: %d%s', status,percent);

        % ------- STATUS UPDATE --------- %
        set(handles.status,'String', str2);
    end

    axes(handles.freq_response);
        semilogx(freq,vectGaindB,'-o','LineWidth',1,'Color',red,'MarkerSize',4);
        ylim ([min(vectGaindB)-marginfreq max(vectGaindB)+marginfreq]); 
        xlim([10 Fs/2*0.9]);
        ylabel('Gain (dB)');
        xlabel('Frequency (Hz)');
        title('Frequency response');
        grid on;


    
% --- Executes on button press in gainbutton.
function gainbutton_Callback(hObject, eventdata, handles)

    %% -------------------------------------------------------------------    
    % RESET
    % -------------------------------------------------------------------- 

    % ---------- STATUS UPDATE -------------- %
    set(handles.status,'String', 'Status: 0%');

    set(handles.gain, 'BackgroundColor', [0.8 0.8 0.8]);
    set(handles.gain, 'String', '');
    amptype     = get(handles.type,'Value'); % 1 = Half-Bridge, 2 = Full-Bridge
    check_sat   = get(handles.checksaturation,'Value');
    THD         = get(handles.THDresult,'Value');

    if  THD == 0.0 && check_sat == 1
        warndlg('Calculate THD to ensure that the signal is not being sent saturated','Check THD');
        return
    elseif THD >= 1.0 && check_sat == 1
        warndlg('THD higher than expected, check the volumes of the PC','Check volumes');
        return
    end


    %% -------------------------------------------------------------------    
    % VARIABLES
    % --------------------------------------------------------------------

    f_in        = str2double(get(handles.fin, 'String'));
    Fs          = str2double(get(handles.editFs,'String'));
    duration    = 1;
    t           = (1/Fs):(1/Fs):duration; 
    Nbits       = 16;
    NChan       = 2;
    A           = 1;
    input       = A*sin(2*pi*f_in*t);


    %% -------------------------------------------------------------------    
    % CALCULATE GAIN
    % --------------------------------------------------------------------

    % Send and receive amplifier's output signal (DA/AD)
    player  = audioplayer(input,Fs);
    data    = audiorecorder(Fs,Nbits,NChan);
    play(player);
    recordblocking(data,duration);
    output  = getaudiodata(data);

    % ---------- STATUS UPDATE --------------- %
    set(handles.status,'String', 'Status: 50%'); 

    out_L       = output(:,1);
    out_R       = output(:,2);

    if amptype == 1 % Half-Bridge
        power_L = out_L'*out_L;
        power_R = out_R'*out_R;
        gain    = power_R/power_L;
        gaindB  = 10*log10(gain);
    else  % Full-Bridge
        output  = out_L - out_R;
        power   = output'*output;
        gain    = power;
        gaindB  = 10*log10(gain);
    end

    str = sprintf('Gain: %.2f dB',gaindB);
    set(handles.gain,'String',str);
    set(handles.gain, 'BackgroundColor', [0.5 0.5 0.5]);
    set(handles.gain, 'ForegroundColor', [1 1 1]);

    % ---------- STATUS UPDATE ---------------- %
    set(handles.status,'String', 'Status: 100%'); 











%% -------------------------------------------------------------------    
% Functions created by default for each element of the interface:
% --------------------------------------------------------------------


function editf_sin_Callback(hObject, eventdata, handles)
% hObject    handle to editf_sin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function editf_sin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editf_sin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editA_sin_Callback(hObject, eventdata, handles)
% hObject    handle to editA_sin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function editA_sin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editA_sin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editFs_Callback(hObject, eventdata, handles)
% hObject    handle to editFs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function editFs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function THDresult_Callback(hObject, eventdata, handles)
function THDresult_CreateFcn(hObject, eventdata, handles)
% hObject    handle to THDresult (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function fs_Callback(hObject, eventdata, handles)
% hObject    handle to fs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fs as text
%        str2double(get(hObject,'String')) returns contents of fs as a double

% --- Executes during object creation, after setting all properties.
function fs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double

% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function numfreq_Callback(hObject, eventdata, handles)
% hObject    handle to numfreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numfreq as text
%        str2double(get(hObject,'String')) returns contents of numfreq as a double

% --- Executes during object creation, after setting all properties.
function numfreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numfreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function fin_Callback(hObject, eventdata, handles)
% hObject    handle to fin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fin as text
%        str2double(get(hObject,'String')) returns contents of fin as a double

% --- Executes during object creation, after setting all properties.
function fin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in checksaturation.
function checksaturation_Callback(hObject, eventdata, handles)
% hObject    handle to checksaturation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checksaturation

% --- Executes on selection change in type.
function type_Callback(hObject, eventdata, handles)
% hObject    handle to type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from type

% --- Executes during object creation, after setting all properties.
function type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function sliderfreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderfreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
