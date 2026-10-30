function udp_sender_close(u)
%UDP_SENDER_CLOSE  Close the UDP sender.
try
    if strcmp(u.kind, 'java')
        u.obj.close();
    else
        delete(u.obj);
    end
catch
end
end
