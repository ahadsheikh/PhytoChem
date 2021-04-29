function _(el) {
    return document.getElementById(el);
}

function uploadFile() {
    const url = _('form-data')['action'];
    const form = _('form-data');

    let formData = new FormData(form);

    let xhr = new XMLHttpRequest();
    xhr.upload.addEventListener("progress", progressHandler, false);
    xhr.addEventListener("load", completeHandler, false);
    xhr.addEventListener("error", errorHandler, false);
    xhr.addEventListener("abort", abortHandler, false);
    xhr.open("POST", url);
    xhr.send(formData);

    // _('cancel-btn').addEventListener('click', e => {
    //     xhr.abort();
    // })
}

function progressHandler(event) {
    _('submit-btn').style.display = 'none';
    _('loading-div').style.display = 'block';
    //_('cancel-btn').style.display = 'block';
    // let percent = (event.loaded / event.total) * 100;
    // console.log(percent);
}

function completeHandler(event) {
    // const message = _('message');
    // message.innerHTML = message.innerHTML + '<p class="text text-danger">Contribution Data Submitted Successfully</p>';

    _('submit-btn').style.display = 'block';
    _('loading-div').style.display = 'none';
    //_('cancel-btn').style.display = 'none';

    // console.log(event.target.responseText);
}

function errorHandler(event) {
    _('message').innerHTML = '<p>Upload Error Happen</p>';
}

function abortHandler(event) {
    _('submit-btn').style.display = 'block';
    _('loading-div').style.display = 'none';
   //  _('cancel-btn').style.display = 'none';

    const message = _('message');
    message.innerHTML = message.innerHTML + '<p class="text text-danger">Upload Aborted</p>';
}

_('form-data').addEventListener('submit', (e) => {
    // e.preventDefault();
    uploadFile();
});
