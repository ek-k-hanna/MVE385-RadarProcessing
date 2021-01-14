function display_time_domain(S, index, PRF)
    IQ = S(:,index)
    time = (0:size(IQ,1) - 1) / PRF;
    plot3(time, real(IQ), imag(IQ))

    %axis equal
    xlabel('Time [s]')
    ylabel('I')
    zlabel('Q')
end