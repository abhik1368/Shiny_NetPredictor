render_helpfile <- function(title, file)
{
    file <- paste0("", file)
    
    body <- markdown::markdownToHTML(file, fragment.only=TRUE, options=c(""))
    link <- paste0(gsub(title, pattern=" ", replacement=""), "_help")
    thisyear <- format(Sys.Date(), "%Y")
    
    html <- sprintf("
<style type='text/css'>
    @media screen and (min-width: 768px) {
        .modal-dialog {
                    width: 800px; /* New width for default modal */
}
</style>
    <hr>
    <div class='modal fade' id='%s' tabindex='-1' role='dialog' aria-labelledby='%s_label' aria-hidden='true' data-width='1000'>
      <div class='modal-dialog'>
        <div class='modal-content'>
          <div class='modal-header'>
            <button type='button' class='close' data-dismiss='modal' aria-label='Close'><span aria-hidden='true'>&times;</span></button>
            <h4 class='modal-title' id='%s_label'>%s</h4>
          </div>
          <div class='modal-body'>%s
          <br><br>
            <font size='1'>
            &copy;%s Abhik Seal.
            All documentation is released under a
              <a rel='license' href='http://creativecommons.org/licenses/by-sa/4.0/' target='_blank'>Creative Commons License</a>.
              <img alt='' style='border-width:0' height=20px src='./img/cc.png'>
            </font>
          </div>
        </div>
      </div>
    </div>
    <div title='Help' data-toggle='modal' data-target='#%s'>
      <i>Click for help</i>
      <div title='Help' class='glyphicon glyphicon-question-sign'></div>
    </div>
     ", link, link, link, title, body, thisyear, link)
    
    html <- enc2utf8(html)
    HTML(html)
}