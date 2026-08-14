using SMTPClient
using CairoMakie


struct MailNotifier

    mailto::Vector{String}
    mailfrom::String
    url::String
    
    send_opts::SendOptions

    function MailNotifier(;mailto::Vector{String}, env_mail_user::String, env_mail_pass::String, url::String)

        mail_user=ENV[env_mail_user]
        mail_pass=ENV[env_mail_pass]
        opt = SendOptions(isSSL = true, username = mail_user, passwd = mail_pass)

        new(mailto, mail_user, url, opt)
    end

    function MailNotifier(;env_mail_user::String, env_mail_pass::String, url::String)

        mail_user=ENV[env_mail_user]
        mail_pass=ENV[env_mail_pass]
        opt = SendOptions(isSSL = true, username = mail_user, passwd = mail_pass)

        mailto = [mail_user]
        new(mailto, mail_user, url, opt)
    end

end

function SMTPClient.send(notifier::MailNotifier, subject::String, message::String)

    
    body = get_body(notifier.mailto, notifier.mailfrom, subject, get_mime_msg(message))

    resp = send(notifier.url, notifier.mailto, notifier.mailfrom, body, notifier.send_opts)

    return resp

end

function SMTPClient.send(notifier::MailNotifier, subject::String, message::String, attachments::Vector{String})

    
    body = get_body(notifier.mailto, notifier.mailfrom, subject, get_mime_msg(message); attachments)

    resp = send(notifier.url, notifier.mailto, notifier.mailfrom, body, notifier.send_opts)

    return resp

end



notifier = MailNotifier(env_mail_user="protonmail_user", env_mail_pass="protonmail_token", url="smtp://smtp.protonmail.ch:587")


session_tmp_dir = mktempdir()


fig, ax, plt = scatter(rand(10), rand(10))


save("$session_tmp_dir/test.pdf", fig)

send(notifier, "Notifier test", "This is a test for the MailNotifier struct, with attachments", ["$session_tmp_dir/test.pdf"])