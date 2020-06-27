function points = loadPerturbedCircle(frequency, amplitude)

    if nargin==0
        amplitude = .025;
        frequency = 100;
    end

    t = linspace(0,2*pi, frequency*20); % enough values of t to guarantee no nyquist sampling issues.
    radiuscenter = amplitude/.05;
    r = radiuscenter + amplitude*sin(frequency*t)';
    
    points = r.*[cos(t)' sin(t)'];
end