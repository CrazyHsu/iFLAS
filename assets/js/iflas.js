$(function() {
  function clickNav(event) {
    var target = event.currentTarget;
    var tab = $(target).attr("href");
    // console.log(tab);
    $('a[href="' + tab + '"]').tab('in');
    $(tab).tab('in');
    if ($("a[role='tab']")) {
      var arrLink = $("a[role='tab']");
      for (i = 0; i < arrLink.length; i++) {
        var link = arrLink[i];
        if (link.classList.contains('active')) {
          link.classList.remove('active');
        }
      }
    }
    target.classList.add('active');
  }

  $(".collapse").on('show.bs.collapse', function() {
    $(this).prev(".sample").find(".icon-plus").removeClass("icon-plus").addClass("icon-minus");
  }).on('hide.bs.collapse', function() {
    $(this).prev(".sample").find(".icon-minus").removeClass("icon-minus").addClass("icon-plus");
    $(this).prev(".sample").removeClass("active");
  });

  if ($("a[role='tab']")) {
    var arrLink = $("a[role='tab']");
    for (i = 0; i < arrLink.length; i++) {
      var link = arrLink[i];
      link.classList.remove('active');
      link.onclick = clickNav;
    }
  }
});
