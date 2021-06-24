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


// Make local date
rows = document.getElementsByClassName("contributer-row")
const months = {
    'January': 1,
    'February': 2,
    'March': 3,
    'April': 4,
    'May': 5,
    'June': 6,
    'July': 7,
    'August': 8,
    'September': 9,
    'October': 10,
    'November': 11,
    'December': 12
}
const monthsRev = ['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December']

for(i = 0; i < rows.length; i++){
    dateTag = rows[i].getElementsByTagName('small')
    UTCTime = dateTag[0].innerText
    let l = UTCTime.split(' ')
    const month = l[0]
    const day = l[1].substring(0, l[1].length-1)
    const year = l[2].substring(0, l[2].length-1)
    const time = l[3]
    const partOfTheDay = l[4].repeat(',','').toUpperCase()
    let lt = new Date(`${month}/${day}/${year} ${time} ${partOfTheDay} UTC`)

    const dateTimeShow = `${lt.getDate()}, ${monthsRev[lt.getMonth()]} ${lt.getFullYear()} ${lt.toLocaleTimeString()}`
    dateTag[0].innerHTML = dateTimeShow
}