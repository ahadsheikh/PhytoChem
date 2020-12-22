// Member dropdown

const $dropdown = $(".dropdown");
const $dropdownToggle = $(".dropdown-toggle");
const $dropdownMenu = $(".dropdown-menu");
const showClass = "show";

$(window).on("load resize", function () {
    if (this.matchMedia("(min-width: 768px)").matches) {
        $dropdown.hover(
            function () {
                const $this = $(this);
                $this.addClass(showClass);
                $this.find($dropdownToggle).attr("aria-expanded", "true");
                $this.find($dropdownMenu).addClass(showClass);
            },
            function () {
                const $this = $(this);
                $this.removeClass(showClass);
                $this.find($dropdownToggle).attr("aria-expanded", "false");
                $this.find($dropdownMenu).removeClass(showClass);
            }
        );
    } else {
        $dropdown.off("mouseenter mouseleave");
    }
});

// Sticky navbar
$(document).ready(function () {
    const pathname = window.location.pathname;
    if (pathname == "/")
        $('.navbar', 'body').addClass("sticky");
    else
        $('.navbar', 'body').removeClass("sticky");
});
$(document).ready(function () {
    const pathname = window.location.pathname;
    if (pathname == "/")
        $('body').addClass("st");
    else
        $('body').removeClass("st");
});
$(document).ready(function () {
    const width = $(window).width();
    const height = $(window).height();
    if (width <= 400 && height <= 800)
        $("#hide").click(function () {
            $(".homepage-text-2").hide();
        });
});

// Change navbar color when scrolling
$(window).scroll(function () {
    const pathname = window.location.pathname;
    if (pathname == "/") {
        if ($(window).scrollTop() >= 50) {
            $('.navbar').css('background', 'linear-gradient(to right,#1a005e, #060009,#330148)');
        } else {
            $('.navbar').css('background', 'transparent');
        }
    }
});